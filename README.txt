REAL-DATA TEST OF compare_a1_frequencies.py
Run date: 2026-11-03

Study
-----
Hof JP and Speed D. LDAK-KVIK performs fast and powerful mixed-model association
analysis of quantitative and binary phenotypes. Nature Genetics 57, 2116-2123
(2025). https://doi.org/10.1038/s41588-025-02286-z

The authors' download page lists UK Biobank standing height, field 50, with
458,079 participants: https://www.ldak-kvik.com/summaries/

Downloaded public study files:
https://genetics.ghpc.au.dk/doug/SumStats/GWAS1.Standing_Height.assoc.gz
https://genetics.ghpc.au.dk/doug/PRS/GWAS1.Standing_Height.effects.gz

Important input detail
----------------------
The downloaded association file contains MAF, but no labelled A1 frequency or
A1 mean dosage. MAF was NOT treated as A1 frequency. The accompanying prediction
model contains Predictor/A1/A2/Centre/Effect. Its Centre is the mean A1 allele
count; for these autosomal SNPs, Centre/2 gives an A1 frequency. This frequency
was aligned to the association file's A1. These are prediction-model frequencies
from the same study's companion file, not a directly reported GWAS EAF column.
The exact equality of the model-training and association-analysis samples was
not independently established. The diagnostic below is explicitly a comparison
using the companion prediction-model frequencies.

The Centre convention is documented at:
https://dougspeed.com/profile-scores/
https://dougspeed.com/blup-estimates/

Reference
---------
1000 Genomes Phase 3, chromosome 1, b37/hs37d5 reference assembly:
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

The VCF's EUR_AF, AFR_AF, EAS_AF, SAS_AF and AMR_AF fields are ALT-allele
frequencies. Reference tables therefore use A1=ALT and A2=REF. Original reference
IDs in the downloaded excerpt are missing. We assigned study rsIDs only after
exact chromosome, position and allele-pair matching. See variant_mapping.tsv;
the original VCF ID is retained there and in the VCF subset. All selected study
coordinates match the b37 reference coordinates. No liftover was performed.

Small test selection
--------------------
Downloaded only bounded compressed prefixes: 1 MiB association output, 1 MiB
prediction output, and 8 MiB reference VCF. All servers returned HTTP 206 with
the requested Content-Range. Headers and byte-range SHA-256 hashes are retained.
The .gz.part files are intentionally incomplete compressed source excerpts;
they are not full-file downloads. The preparer reads complete records and stops
at the expected truncated end. The saved .subset files are complete small files.

Selected the first 200 PASS biallelic SNPs in reference VCF order with overlap
in both study files and matching positions/alleles. No frequency-similarity
criterion was used to choose SNPs. The selected interval is chromosome 1,
752721-1555718 (b37). All five reference tables contain exactly the same SNPs.
The comparison script excludes 8 A/T or C/G variants, leaving 192 SNPs.

Results
-------
Panel  Matched SNPs  Pearson r   Mean absolute difference  RMSE
EUR    192           0.999423    0.008081                  0.010869
AMR    192           0.944400    0.065457                  0.106105
SAS    192           0.919227    0.078949                  0.130215
EAS    192           0.759089    0.138102                  0.233271
AFR    192           0.627850    0.171550                  0.267282

Frequency differences use the 0-1 scale; 0.008081 is 0.8081 percentage points.
EUR is the closest of these five references in this small subset. The variants
are from one genomic interval and are not independent genome-wide markers.
These results demonstrate the calculation on real study data; they do not
establish a general ancestry-matching cutoff or validate LDSC performance.

Validation
----------
The actual standalone script was invoked separately for all five panels.
Every result agreed with an independent two-pass calculation using Python's
statistics.correlation and math.fsum (relative tolerance 1e-12, absolute
tolerance 1e-14). Centre/2 and allele orientation were checked against the
original prediction records. Full-precision results, script SHA-256 and the
Python version are in results.json. The script required no modification.

Rerun using the included prepared data (no network needed)
--------------------------------------------------------
python3 "/Users/breezece/Google Drive official async/ldscore.ms/ldak_frequency_example/run_example.py"

Run the minimal comparator for EUR only:
python3 "/Users/breezece/Google Drive official async/ldscore.ms/compare_a1_frequencies.py" \
  "/Users/breezece/Google Drive official async/ldscore.ms/ldak_frequency_example/height.frequencies.tsv" \
  "/Users/breezece/Google Drive official async/ldscore.ms/ldak_frequency_example/1000g.EUR.frequencies.tsv"

Rebuild the small tables from the retained downloaded prefixes:
python3 "/Users/breezece/Google Drive official async/ldscore.ms/ldak_frequency_example/prepare_example.py"

Optional fresh download (overwrites the three retained prefixes and headers):
bash "/Users/breezece/Google Drive official async/ldscore.ms/ldak_frequency_example/download_example.sh"

Python 3.10+ is sufficient for the example runner. No external Python packages
are required. File names and exact preprocessing are recorded in the scripts;
results.tsv is a compact tabular copy of the results.
