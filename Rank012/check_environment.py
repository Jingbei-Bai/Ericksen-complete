"""Check the exact Linux environment used for the independent recomputation."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import platform
import sys
from urllib.parse import urlsplit, unquote


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error('Use Sage Python with assertions enabled')
    here = Path(__file__).resolve().parent
    lock = (here/'environment-linux-64.lock').read_text()
    expected = {unquote(Path(urlsplit(line).path).name): urlsplit(line).fragment
                for line in lock.splitlines() if line.startswith('https://')}
    metadata = Path(sys.prefix)/'conda-meta'
    actual = {}
    for path in metadata.glob('*.json'):
        p = json.loads(path.read_text())
        actual[p['fn']] = p['md5']
    differences = {name: {'expected_md5': value, 'actual_md5': actual.get(name)}
                   for name, value in expected.items() if actual.get(name) != value}
    unexpected_packages=sorted(set(actual)-set(expected))
    import sage.version
    import sympy
    from sage.all import PolynomialRing, QQ
    versions = dict(sage=sage.version.version, python=platform.python_version(),
                    sympy=sympy.__version__)
    required = dict(sage='10.9', python='3.12.14', sympy='1.14.0')
    assert versions == required, (versions, required)
    assert platform.system() == 'Linux' and platform.machine() == 'x86_64'
    assert not differences, differences
    assert not unexpected_packages, unexpected_packages
    # Exercise the exact multivariate ideal arithmetic used by the certificates.
    ring = PolynomialRing(QQ, names=('x', 'y'))
    x, y = ring.gens()
    ideal = ring.ideal([x*y-1, x**2-y])
    assert ideal.reduce(x**3-1) == 0
    fingerprint = hashlib.sha256(json.dumps(sorted(actual.items())).encode()).hexdigest()
    record = dict(schema='ericksen_environment_check/v1', status='PASS',
                  checked_utc=datetime.now(timezone.utc).isoformat(),
                  versions=versions, platform='linux-64',
                  conda_package_count=len(actual), required_package_count=len(expected),
                  mismatches=differences, package_fingerprint_sha256=fingerprint,
                  unexpected_packages=unexpected_packages,
                  executable=sys.executable, prefix=sys.prefix,
                  lock_sha256=hashlib.sha256((here/'environment-linux-64.lock').read_bytes()).hexdigest(),
                  exact_ideal_arithmetic_check=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2), encoding='utf-8')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
