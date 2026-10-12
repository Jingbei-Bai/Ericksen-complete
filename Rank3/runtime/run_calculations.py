"""Replay an explicit calculation DAG with the active Python interpreter.

Use SageMath's Python for execution. --list only reads and displays the DAG.
The companion execution_jobs.json contains the mathematical job definitions.
"""
from __future__ import annotations

import argparse
import heapq
import json
import os
from pathlib import Path, PureWindowsPath
import shutil
import subprocess
import sys


class JobError(Exception):
    """An invalid job definition or a failed replay action."""


def member_path(root: Path, member: str) -> Path:
    if not isinstance(member, str) or not member:
        raise JobError("A bundle member must be a nonempty string.")
    normalized = member.replace("\\", "/")
    if Path(normalized).is_absolute() or PureWindowsPath(normalized).drive:
        raise JobError(f"Bundle member must be relative: {member!r}")
    path = (root / normalized).resolve()
    if not path.is_relative_to(root):
        raise JobError(f"Bundle member escapes the bundle root: {member!r}")
    return path


def read_jobs(path: Path) -> list[dict]:
    try:
        data = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, ValueError) as exc:
        raise JobError(f"Cannot read {path}: {exc}") from exc
    if not isinstance(data, dict) or not isinstance(data.get("jobs"), list):
        raise JobError("execution_jobs.json must contain a jobs array.")
    jobs = data["jobs"]
    seen = set()
    for number, job in enumerate(jobs, 1):
        if not isinstance(job, dict):
            raise JobError(f"Job {number} must be an object.")
        job_id = job.get("id")
        if not isinstance(job_id, str) or not job_id.strip():
            raise JobError(f"Job {number} requires a nonempty string id.")
        if job_id in seen:
            raise JobError(f"Duplicate job id: {job_id}")
        seen.add(job_id)
        if not isinstance(job.get("group"), str) or not job["group"]:
            raise JobError(f"Job {job_id} requires a nonempty group.")
        deps = job.get("depends_on", [])
        if not isinstance(deps, list) or any(not isinstance(x, str) for x in deps):
            raise JobError(f"Job {job_id}: depends_on must be a string array.")
        if len(deps) != len(set(deps)):
            raise JobError(f"Job {job_id} repeats a dependency.")
        args = job.get("arguments", [])
        if not isinstance(args, list) or any(not isinstance(x, str) for x in args):
            raise JobError(f"Job {job_id}: arguments must be a string array.")
        action = job.get("action", "run")
        if action == "run":
            if not isinstance(job.get("source_member"), str) or not job["source_member"]:
                raise JobError(f"Run job {job_id} requires source_member.")
            wrapper = job.get("wrapper_member")
            if wrapper is not None and (not isinstance(wrapper, str) or not wrapper):
                raise JobError(f"Job {job_id}: wrapper_member must be a nonempty string or null.")
        elif action == "copy_generated":
            if any(not isinstance(job.get(key), str) or not job[key] for key in ("from", "to")):
                raise JobError(f"Copy job {job_id} requires from and to members.")
        else:
            raise JobError(f"Job {job_id} has unsupported action {action!r}.")
    return jobs


def topological_jobs(jobs: list[dict], groups: list[str], job_ids: list[str] | None = None) -> list[dict]:
    by_id = {job["id"]: job for job in jobs}
    positions = {job["id"]: i for i, job in enumerate(jobs)}
    following = {job_id: [] for job_id in by_id}
    degrees = {}
    for job in jobs:
        deps = job.get("depends_on", [])
        missing = [dep for dep in deps if dep not in by_id]
        if missing:
            raise JobError(f"Job {job['id']} has unknown dependencies: {', '.join(missing)}")
        degrees[job["id"]] = len(deps)
        for dep in deps:
            following[dep].append(job["id"])
    ready = [(positions[job_id], job_id) for job_id, n in degrees.items() if n == 0]
    heapq.heapify(ready)
    ordered = []
    while ready:
        _, job_id = heapq.heappop(ready)
        ordered.append(by_id[job_id])
        for child in following[job_id]:
            degrees[child] -= 1
            if degrees[child] == 0:
                heapq.heappush(ready, (positions[child], child))
    if len(ordered) != len(jobs):
        blocked = [job_id for job_id, n in degrees.items() if n]
        raise JobError("Dependency cycle blocks jobs: " + ", ".join(blocked))
    if not groups and not job_ids:
        return ordered
    unknown = set(groups) - {job["group"] for job in jobs}
    if unknown:
        raise JobError("Unknown groups: " + ", ".join(sorted(unknown)))
    unknown_jobs = set(job_ids or []) - set(by_id)
    if unknown_jobs:
        raise JobError("Unknown jobs: " + ", ".join(sorted(unknown_jobs)))
    selected = {job["id"] for job in jobs if job["group"] in groups} | set(job_ids or [])
    pending = list(selected)
    while pending:
        for dep in by_id[pending.pop()].get("depends_on", []):
            if dep not in selected:
                selected.add(dep)
                pending.append(dep)
    return [job for job in ordered if job["id"] in selected]


def command(job: dict, root: Path) -> list[str]:
    executable = member_path(root, job.get("wrapper_member") or job["source_member"])
    args = []
    for value in job.get("arguments", []):
        normalized = value.replace("\\", "/")
        args.append(str(member_path(root, normalized)) if normalized.startswith(("research/", "output/")) else value)
    runtime = member_path(root, "runtime_execution.py")
    return [sys.executable, str(runtime), str(executable), *args]


def validate_runtime(jobs: list[dict], root: Path) -> None:
    for job in jobs:
        if job.get("action", "run") == "copy_generated":
            member_path(root, job["from"])
            member_path(root, job["to"])
            continue  # A preceding producer may create the copy source.
        members = [job["source_member"]]
        if job.get("wrapper_member"):
            members.append(job["wrapper_member"])
        for member in members:
            source = member_path(root, member)
            if not source.is_file():
                raise JobError(f"Job {job['id']}: runtime source is missing: {source}")
        command(job, root)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="Display ordered jobs without executing or copying files.")
    parser.add_argument("--group", action="append", default=[], metavar="GROUP", help="Select a group and all its prerequisites; repeat for multiple groups.")
    parser.add_argument("--job", action="append", default=[], metavar="ID", help="Select a job and all its prerequisites; repeat for multiple jobs.")
    options = parser.parse_args(argv)
    script_dir = Path(__file__).resolve().parent
    root = Path(os.environ.get("ERICKSEN_BUNDLE_ROOT", str(script_dir))).expanduser().resolve()
    try:
        jobs = topological_jobs(read_jobs(script_dir / "execution_jobs.json"), options.group, options.job)
        if options.list:
            for number, job in enumerate(jobs, 1):
                action = job.get("action", "run")
                detail = f"{job['from']} -> {job['to']}" if action == "copy_generated" else json.dumps(command(job, root), ensure_ascii=False)
                print(f"{number:03d} {job['id']} [{job['group']}] {action}: {detail}")
            return 0
        validate_runtime(jobs, root)
        environment = dict(os.environ, ERICKSEN_BUNDLE_ROOT=str(root))
        for number, job in enumerate(jobs, 1):
            print(f"START {number}/{len(jobs)} {job['id']} [{job['group']}]", flush=True)
            if job.get("action", "run") == "copy_generated":
                source, target = member_path(root, job["from"]), member_path(root, job["to"])
                if not source.is_file():
                    raise JobError(f"Job {job['id']}: generated copy source is missing: {source}")
                if source == target:
                    raise JobError(f"Job {job['id']}: copy source and destination are identical.")
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, target)
            else:
                completed = subprocess.run(command(job, root), cwd=str(root), env=environment)
                if completed.returncode:
                    print(f"FAIL {job['id']}: child exit code {completed.returncode}", file=sys.stderr, flush=True)
                    return completed.returncode if 0 < completed.returncode < 256 else 1
            print(f"DONE {job['id']}", flush=True)
        print(f"Completed {len(jobs)} calculation jobs.", flush=True)
        return 0
    except (JobError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
