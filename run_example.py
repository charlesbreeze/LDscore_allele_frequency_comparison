#!/usr/bin/env python3
"""Run the standalone script and independently verify its real-data results."""
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT.parent / 'compare_a1_frequencies.py'
PANELS = ('EUR', 'AFR', 'EAS', 'SAS', 'AMR')


def read_table(name):
    with (ROOT / name).open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def main():
    gwas = read_table('height.frequencies.tsv')
    effects = {r['Predictor']: r for r in read_table('height.effects.subset.tsv')}
    # Check the conversion against the original prediction-model records.
    for row in gwas:
        model = effects[row['SNP']]
        frequency = float(model['Centre']) / 2
        if row['A1'] != model['A1']:
            assert row['A1'] == model['A2'] and row['A2'] == model['A1']
            frequency = 1 - frequency
        assert math.isclose(float(row['EAF']), frequency, abs_tol=1e-15)

    results = {}
    for panel in PANELS:
        command = [sys.executable, str(SCRIPT), str(ROOT / 'height.frequencies.tsv'),
                   str(ROOT / f'1000g.{panel}.frequencies.tsv')]
        completed = subprocess.run(command, capture_output=True, text=True, check=True)
        result = json.loads(completed.stdout)
        reference = {r['SNP']: r for r in read_table(f'1000g.{panel}.frequencies.tsv')}
        assert set(reference) == {r['SNP'] for r in gwas}
        pairs = []
        for row in gwas:
            if {row['A1'], row['A2']} in ({'A', 'T'}, {'C', 'G'}):
                continue
            ref = reference[row['SNP']]
            assert {row['A1'], row['A2']} == {ref['A1'], ref['A2']}
            frequency = float(ref['EAF'])
            if row['A1'] != ref['A1']:
                frequency = 1 - frequency
            pairs.append((float(row['EAF']), frequency))
        x, y = zip(*pairs)
        expected = {
            'matched_snps': len(pairs),
            'pearson_r': statistics.correlation(x, y),
            'mean_absolute_difference': math.fsum(abs(a - b) for a, b in pairs) / len(pairs),
            'rmse': math.sqrt(math.fsum((a - b) ** 2 for a, b in pairs) / len(pairs))
        }
        for key, value in expected.items():
            if not math.isclose(result[key], value, rel_tol=1e-12, abs_tol=1e-14):
                raise AssertionError(f'{panel}: disagreement in {key}')
        results[panel] = result
    assert len({r['matched_snps'] for r in results.values()}) == 1
    report = {
        'results': results,
        'validation': 'All five CLI outputs agree with independent two-pass statistics; Centre/2 conversion checked against original prediction records',
        'script_sha256': hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        'python_version': sys.version,
    }
    (ROOT / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    with (ROOT / 'results.tsv').open('w') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        columns = ['matched_snps', 'pearson_r', 'mean_absolute_difference', 'rmse']
        writer.writerow(['reference_panel', *columns])
        for panel, result in sorted(results.items(), key=lambda pair: pair[1]['mean_absolute_difference']):
            writer.writerow([panel, *(result[c] for c in columns)])
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
