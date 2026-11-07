#!/usr/bin/env bash
# Download bounded prefixes, retaining the server's Content-Range and version metadata.
set -euo pipefail
cd -- "$(dirname -- "$0")"
mkdir -p downloads
curl -L --fail --silent --show-error --max-time 120 \
  --range 0-1048575 --max-filesize 1048576 \
  -D downloads/height.headers.txt \
  'https://genetics.ghpc.au.dk/doug/SumStats/GWAS1.Standing_Height.assoc.gz' \
  -o downloads/height.assoc.gz.part
curl -L --fail --silent --show-error --max-time 120 \
  --range 0-1048575 --max-filesize 1048576 \
  -D downloads/height_effects.headers.txt \
  'https://genetics.ghpc.au.dk/doug/PRS/GWAS1.Standing_Height.effects.gz' \
  -o downloads/height.effects.gz.part
curl -L --fail --silent --show-error --max-time 120 \
  --range 0-8388607 --max-filesize 8388608 \
  -D downloads/1000g.headers.txt \
  'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz' \
  -o downloads/1000g.chr1.vcf.gz.part
