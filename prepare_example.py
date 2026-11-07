#!/usr/bin/env python3
"""Prepare a small reproducible comparison from downloaded public file prefixes.

Run after download_example.sh. Only complete newline-terminated records are used.
The published association file lacks A1 frequency; frequency is obtained from
the accompanying LDAK prediction model's Centre/2, then aligned to GWAS A1.
"""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DOWNLOADS = ROOT / 'downloads'
PANELS = ('EUR', 'AFR', 'EAS', 'SAS', 'AMR')


def complete_lines(path):
    # These are intentional HTTP byte-range prefixes, not complete gzip files.
    with gzip.open(path, 'rt') as stream:
        while True:
            try:
                line = next(stream)
            except (StopIteration, EOFError):
                return
            if not line.endswith('\n'):
                return
            yield line.rstrip('\n')


def read_table(name):
    lines = complete_lines(DOWNLOADS / name)
    header = next(lines).split()
    records = {}
    for line in lines:
        values = line.split()
        if len(values) != len(header):
            raise ValueError(f'Malformed record in {name}')
        row = dict(zip(header, values))
        snp = row['Predictor']
        if snp in records:
            raise ValueError(f'Duplicate SNP in {name}: {snp}')
        records[snp] = row
    return header, records


def write_table(name, header, rows):
    with (ROOT / name).open('w') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


def main():
    assoc_header, assoc = read_table('height.assoc.gz.part')
    effects_header, effects = read_table('height.effects.gz.part')
    by_position = {}
    for snp, row in assoc.items():
        if snp not in effects:
            continue
        key = (row['Chromosome'], row['Basepair'], frozenset((row['A1'], row['A2'])))
        if key in by_position:
            raise ValueError(f'Ambiguous study coordinate/allele key: {key}')
        by_position[key] = snp
    selected, sites, vcf_meta = [], [], []
    seen = set()
    for line in complete_lines(DOWNLOADS / '1000g.chr1.vcf.gz.part'):
        if line.startswith('##'):
            vcf_meta.append(line)
            continue
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        chrom, pos, original_id, ref, alt, qual, status, info = fields[:8]
        if status != 'PASS' or ref not in 'ACGT' or alt not in 'ACGT':
            continue
        if len(ref) != 1 or len(alt) != 1 or ref == alt:
            continue
        # This reference release has missing IDs; use exact coordinates AND alleles.
        snp = by_position.get((chrom, pos, frozenset((ref, alt))))
        if snp is None or snp in seen:
            continue
        a, e = assoc[snp], effects[snp]
        if a['Chromosome'] != chrom or a['Basepair'] != pos:
            raise ValueError(f'Coordinate mismatch for {snp}')
        if {a['A1'], a['A2']} != {ref, alt}:
            continue
        frequency = float(e['Centre']) / 2
        if not math.isfinite(frequency) or not 0 <= frequency <= 1:
            raise ValueError(f'Invalid Centre for {snp}')
        if (e['A1'], e['A2']) == (a['A1'], a['A2']):
            pass
        elif (e['A2'], e['A1']) == (a['A1'], a['A2']):
            frequency = 1 - frequency
        else:
            raise ValueError(f'Prediction/GWAS allele mismatch for {snp}')
        attributes = dict(item.split('=', 1) for item in info.split(';') if '=' in item)
        frequencies = {p: float(attributes[p + '_AF']) for p in PANELS}
        if any(not math.isfinite(f) or not 0 <= f <= 1 for f in frequencies.values()):
            raise ValueError(f'Invalid reference frequency for {snp}')
        selected.append((snp, frequency, alt, ref, frequencies))
        sites.append('\t'.join(fields[:8]))
        seen.add(snp)
        if len(selected) == 200:
            break
    if len(selected) != 200:
        raise ValueError(f'Need 200 overlapping SNPs; found {len(selected)}')

    write_table('height.assoc.subset.tsv', assoc_header,
                ([assoc[s][h] for h in assoc_header] for s, *_ in selected))
    write_table('height.effects.subset.tsv', effects_header,
                ([effects[s][h] for h in effects_header] for s, *_ in selected))
    write_table('height.frequencies.tsv', ['SNP', 'A1', 'A2', 'EAF'],
                ([s, assoc[s]['A1'], assoc[s]['A2'], format(f, '.17g')]
                 for s, f, *_ in selected))
    write_table('variant_mapping.tsv', ['SNP', 'CHR', 'BP', 'REF', 'ALT', 'ORIGINAL_VCF_ID'],
                ([s, assoc[s]['Chromosome'], assoc[s]['Basepair'], ref, alt,
                  sites[i].split('\t')[2]]
                 for i, (s, _, alt, ref, _) in enumerate(selected)))
    for panel in PANELS:
        write_table(f'1000g.{panel}.frequencies.tsv', ['SNP', 'A1', 'A2', 'EAF'],
                    ([s, alt, ref, format(freqs[panel], '.17g')]
                     for s, _, alt, ref, freqs in selected))
    with gzip.open(ROOT / '1000g.chr1.sites.subset.vcf.gz', 'wt') as stream:
        stream.write('\n'.join(vcf_meta) + '\n')
        stream.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        stream.write('\n'.join(sites) + '\n')
    metadata = {
        'selected_snps': len(selected),
        'selection': 'First 200 PASS biallelic SNPs in VCF order present in both study files, with matching chromosome, position and allele pair; no frequency-similarity selection',
        'chromosome': '1',
        'first_position': int(assoc[selected[0][0]]['Basepair']),
        'last_position': int(assoc[selected[-1][0]]['Basepair']),
        'coordinate_check': 'All selected study positions match the b37 reference VCF positions',
        'study_frequency_source': 'LDAK prediction model Centre/2, aligned to association A1; not MAF from the association file',
        'reference_frequency_source': 'VCF EUR_AF/AFR_AF/EAS_AF/SAS_AF/AMR_AF fields; reference table A1=ALT, A2=REF',
        'reference_identifier_mapping': 'Original VCF IDs are missing; assign study rsIDs only after exact chromosome, position and allele-pair matching; see variant_mapping.tsv',
        'input_prefix_complete_records': {'association': len(assoc), 'prediction': len(effects)},
        'sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in sorted(DOWNLOADS.glob('*.part'))}
    }
    (ROOT / 'preparation.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
