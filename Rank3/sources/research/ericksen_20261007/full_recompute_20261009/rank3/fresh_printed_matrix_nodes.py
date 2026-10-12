"""Read the current supplement and compare its 20 moment entries with fresh sources.

Run with SageMath's Python in the extracted proof bundle. --inspect-tex uses
ordinary Python to check the TeX declarations without performing arithmetic.
The source derivation must first generate fresh_source_to_moments.sobj.
Row normalization factors are read from the current supplement itself.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import time

RESEARCH_REL = Path('research/ericksen_20261007')


def read_tex_tree(path, paper_root, active=()):
    """Expand only local input files; keep their paths for the source hashes."""
    path = path.resolve()
    if not path.is_relative_to(paper_root.resolve()):
        raise ValueError('TeX input escapes the paper directory: ' + str(path))
    if path in active:
        raise ValueError('Recursive TeX input: ' + str(path))
    raw = path.read_text(encoding='utf-8')
    records = {path: hashlib.sha256(path.read_bytes()).hexdigest()}

    def expand(match):
        target = paper_root / match[1]
        if not target.suffix:
            target = target.with_suffix('.tex')
        text, hashes = read_tex_tree(target, paper_root, (*active, path))
        records.update(hashes)
        return text

    return re.sub(r'\\input\{([^}]+)\}', expand, raw), records


def unwrap(text):
    return text.replace(r'\begin{aligned}', '').replace(r'\end{aligned}', '')


def fractions(text):
    def group(start):
        if start >= len(text) or text[start] != '{':
            raise ValueError('Expected braced fraction argument')
        depth = 1
        end = start + 1
        while end < len(text) and depth:
            depth += (text[end] == '{') - (text[end] == '}')
            end += 1
        if depth:
            raise ValueError('Unclosed fraction argument')
        return text[start + 1:end - 1], end

    while r'\frac' in text:
        start = text.index(r'\frac')
        numerator, end = group(start + len(r'\frac'))
        denominator, end = group(end)
        value = '(' + fractions(numerator) + ')/(' + fractions(denominator) + ')'
        text = text[:start] + value + text[end:]
    return text


def arithmetic(text):
    """The original polynomial-only tokenizer, with current symbol aliases."""
    text = unwrap(text).replace(r'\left', '').replace(r'\right', '')
    text = text.replace(r'\qquad', ' ').replace(r'\quad', ' ')
    text = text.replace('\\\\', '').replace('&', '')
    for macro, name in ((r'\varsigma', 'S'), (r'\eta', 'H'), (r'\zeta', 'Z')):
        text = text.replace(macro, ' ' + name + ' ')
    text = re.sub(r'\\gamma_\{(\d+)\}', lambda m: ' g' + m[1] + ' ', text)
    text = fractions(text)
    text = re.sub(r'\^\{(\d+)\}', lambda m: '**' + m[1], text).replace('^', '**')
    if re.search(r'[\\{}]', text):
        raise ValueError('Unsupported polynomial TeX syntax: ' + text[:100])
    tokens = re.findall(r'\*\*|[A-Za-z][A-Za-z0-9_]*|\d+|[()+\-*/]', text)
    if ''.join(tokens) != re.sub(r'\s+', '', text):
        raise ValueError('Incomplete polynomial tokenization: ' + repr(text[:150]))
    converted = []
    atom = r'[A-Za-z][A-Za-z0-9_]*|\d+'
    for token in tokens:
        if converted and (re.fullmatch(atom, converted[-1]) or converted[-1] == ')') \
                and (re.fullmatch(atom, token) or token == '('):
            converted.append('*')
        converted.append(token)
    return ''.join(converted)


def declarations(tex):
    nodes, entries, row_units = [], {}, {}
    for block in re.findall(r'\\\[(.*?)\\\]', tex, flags=re.S):
        clean = unwrap(block).lstrip()
        if clean.startswith(r'\gamma_'):
            matches = list(re.finditer(r'\\gamma_\{(\d+)\}\s*&?=', clean))
            for position, match in enumerate(matches):
                if int(match[1]) != len(nodes):
                    raise ValueError('Non-topological printed gamma declaration')
                end = matches[position+1].start() if position+1 < len(matches) else len(clean)
                expression = arithmetic(clean[match.end():end])
                dependencies = [int(n) for n in re.findall(r'\bg(\d+)\b', expression)]
                if any(n >= len(nodes) for n in dependencies):
                    raise ValueError('Gamma uses an undefined node')
                nodes.append(expression)
        elif clean.startswith(r'\mathsf A_'):
            match = re.fullmatch(r'\s*\\mathsf A_\{(\d+),(\d+)\}\s*&?=(.*)', clean, flags=re.S)
            if match is None:
                raise ValueError('Invalid printed matrix declaration')
            key = (int(match[1]) - 1, int(match[2]) - 1)
            if key in entries:
                raise ValueError('Duplicate printed matrix entry')
            entries[key] = arithmetic(match[3])
        elif r'\mathrm{row}' in clean and r'd_{' in clean:
            clean = re.sub(r'd_\{([1-5])\}\^\{\\mathrm\{row\}\}', r'D\1', clean)
            if not clean.lstrip().startswith('D'):
                continue
            for piece in re.split(r'\\\\|\\qquad', clean):
                piece = re.sub(r'[\s&]', '', piece).rstrip(',.')
                match = re.fullmatch(r'(D[1-5](?:=D[1-5])*)=(.*)', piece, flags=re.S)
                if match is None:
                    raise ValueError('Invalid printed row factor')
                expression = arithmetic(match[2])
                for name in match[1].split('='):
                    row_units[int(name[1]) - 1] = expression
    if len(nodes) != 425 or set(entries) != {(i, j) for i in range(5) for j in range(4)}:
        raise ValueError('Expected 425 gamma nodes and all 20 matrix entries')
    if set(row_units) != set(range(5)):
        raise ValueError('Expected all five printed row factors')
    return nodes, entries, row_units


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle-root', type=Path)
    parser.add_argument('--paper-root', type=Path)
    parser.add_argument('--inspect-tex', action='store_true')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    bundle = args.bundle_root.resolve() if args.bundle_root is not None else None
    if bundle is None and args.paper_root is None:
        bundle = next((p for p in Path(__file__).resolve().parents
                       if (p / 'paper/supplement.tex').is_file()), None)
        if bundle is None:
            parser.error('Run inside the extracted proof bundle or pass --bundle-root')
    paper = (args.paper_root or bundle / 'paper').resolve()
    entry = paper / 'supplement.tex'
    if not entry.is_file():
        entry = paper / 'Computational_supplement.tex'
    tex, hashes = read_tex_tree(entry, paper)
    nodes, entries, printed_units = declarations(tex)
    if args.inspect_tex:
        print(json.dumps(dict(status='PARSE_PASS',gamma_nodes=len(nodes),matrix_entries=len(entries),
                              printed_row_factors=len(printed_units),tex_files=len(hashes),
                              arithmetic_executed=False), indent=2))
        return

    from sage.all import PolynomialRing, QQ, load, sage_eval
    start = time.monotonic()
    checks = []
    research = bundle / RESEARCH_REL
    out = research / 'full_recompute_20261009/rank3'
    C = PolynomialRing(QQ, names=('H', 'Z', 'S', 'w'))
    H, Z, S, w = C.gens()
    K = C.fraction_field()
    env = dict(H=H, Z=Z, S=S, w=w, L=H*S-Z)

    def check(name, value):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    for i, expression in enumerate(nodes):
        env['g' + str(i)] = sage_eval(expression, locals=env, preparse=False)
    values = {key: sage_eval(expression, locals=env, preparse=False)
              for key, expression in entries.items()}
    fresh = load(str(out / 'fresh_source_to_moments.sobj'))['source_matrix']
    oldK = fresh.base_ring().base_ring()
    phi = oldK.ring().hom([K(H), K(Z), K(S)], K)

    def coefficient(c):
        return phi(c.numerator()) / phi(c.denominator())

    for i in range(5):
        unit = sage_eval(printed_units[i], locals=env, preparse=False)
        check('printed row factor is nonzero ' + str(i), K(unit) != 0)
        for j in range(4):
            target = sum(coefficient(c)*w**int(k) for k, c in fresh[i, j].dict().items())
            check('printed expression actual source ' + str((i, j)),
                  K(values[i, j])*unit == target)
    report = dict(status='PASS', checks=checks, count=len(checks), printed_nodes=len(nodes),
                  printed_entries=len(entries), printed_row_factors=len(printed_units),
                  seconds=time.monotonic()-start, specializations=[],
                  source='Current supplement TeX compared with the freshly source-derived moment matrix',
                  tex_sha256={p.relative_to(paper).as_posix():h for p, h in hashes.items()})
    report_path = args.report or out / 'fresh_printed_matrix_nodes.json'
    report_path.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
