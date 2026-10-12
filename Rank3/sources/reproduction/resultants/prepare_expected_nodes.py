"""Prepare exact candidate values for the unchanged full resultant reader.

These files contain evaluations of the supplied complete candidate polynomials.
They are not independent resultant certificates. The full original reader must
separately recompute every ZZ resultant, compare it to these evaluations, and
verify its Sylvester degree bound and enough distinct exact evaluation points.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import sys
import time
from sage.all import ZZ, PolynomialRing, load

NAMES = (
    'independent_ratio_canonical_second_resultant12',
    'independent_ratio_canonical_second_resultant13',
    'independent_ratio_canonical_second_extra_resultant_index31',
)
NODEFILES = (NAMES[0] + '_interpolation.nodes',
             NAMES[1] + '_interpolation.nodes', NAMES[2] + '.nodes')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    options = parser.parse_args()
    sys.set_int_max_str_digits(0)
    root = options.input_root.resolve()
    out = options.output_root.resolve()
    out.mkdir(parents=True, exist_ok=True)
    P = PolynomialRing(ZZ, 'expected_r')
    records = []
    started = time.monotonic()
    for name, filename in zip(NAMES, NODEFILES):
        source = root / (name + '.sobj')
        payload = load(str(source))
        candidate = P(payload['resultant'])
        assert candidate.degree() <= 2652
        assert candidate != 0
        target = out / filename
        temporary = target.with_suffix(target.suffix + '.tmp')
        if target.exists() or temporary.exists():
            raise FileExistsError('Use an empty generated node destination: ' + str(target))
        digest = hashlib.sha256()
        with temporary.open('w', encoding='ascii', newline='\n') as handle:
            for node in range(2, 2655):
                actual_node = ZZ(node)
                expected = ZZ(candidate(actual_node))
                line = str(actual_node) + ' ' + str(expected) + '\n'
                handle.write(line)
                digest.update(line.encode('ascii'))
        os.replace(temporary, target)
        assert sha(target) == digest.hexdigest()
        records.append(dict(
            candidate_input=str(source), candidate_sha256=sha(source),
            candidate_degree=int(candidate.degree()), generated_file=str(target),
            sha256=digest.hexdigest(), bytes=target.stat().st_size,
            nodes=2653, first_node=2, last_node=2654,
            calculation='Exact complete ZZ-polynomial candidate evaluation only',
        ))
    report = dict(
        schema='ericksen_exact_resultant_expected_node_preparation/v1',
        status='PREPARED_EXACT_CANDIDATE_VALUES_NOT_RESULTANT_VERIFICATION',
        old_node_files_read=False, output_records=records,
        seconds=time.monotonic() - started,
        required_complete_verification='Run the unchanged independent_canonical_main_resultant_reader.py: all7959 full exact resultants, candidate equality, rigorous2652 degree bounds, protected stripping and exact gcd.',
    )
    (out / 'canonical_resultant_expected_node_preparation.json').write_text(
        json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(report['status'], 'all7959 candidate evaluations prepared', flush=True)

if __name__ == '__main__':
    main()
