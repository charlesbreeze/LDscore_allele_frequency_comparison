# Comparing allele frequencies with a reference panel

`compare_a1_frequencies.py` is a Python script that compares how common genetic variants are in a study and in a reference population. It can help researchers assess whether a reference panel is a reasonable candidate for their analysis and spot possible problems with allele labels.

A **reference panel** is a collection of genetic data used for comparison. A **SNP** is a position in the genome where people can have different DNA letters, called **alleles**. An **allele frequency** tells us how common a particular allele is. For example, a frequency of `0.20` means that 20% of the observed allele copies carry that allele.

The script matches SNPs between two files, makes sure their frequencies refer to the same allele, and reports four measures of similarity. It runs locally using Python's standard library, with no additional packages or internet connection needed once the input files are available.

## Try the included example

Use Python 3.10 or newer. In a terminal, run:

```bash
python3 "compare_a1_frequencies.py" \
  "height.frequencies.tsv" \
  "1000g.EUR.frequencies.tsv"
```

This compares frequencies derived from a published LDAK-KVIK height study's companion prediction model with the European group in the 1000 Genomes reference panel. The files are already included in this folder. If you move the folder, update the paths in the command.

The result is printed in the terminal as JSON, a labelled text format. Rounded for readability, it is:

```json
{
  "matched_snps": 192,
  "pearson_r": 0.999423,
  "mean_absolute_difference": 0.008081,
  "rmse": 0.010869
}
```

## Understanding the results

| Output | Meaning | How to read it |
|---|---|---|
| `matched_snps` | Number of SNPs actually compared after matching and filtering. | Check this before interpreting the other results. A small or unrepresentative set may give a misleading impression of overall similarity. |
| `pearson_r` | Correlation between the study and reference allele frequencies. | Values closer to 1 indicate that the frequencies follow a similar pattern. Correlation alone does not measure how close their actual values are. |
| `mean_absolute_difference` | Average absolute difference between the two frequencies. | Smaller is closer. A value of `0.01` means an average difference of **one percentage point**. |
| `rmse` | Root mean square error: another measure of frequency differences that gives more weight to large discrepancies. | Smaller is closer. It uses the same 0–1 frequency scale. |

For example, frequencies of `0.20` and `0.23` differ by `0.03`, or three percentage points. In the included European comparison, the average absolute difference is about **0.81 percentage points**.

If there are fewer than two matched SNPs, or all retained frequencies in either file are identical, correlation cannot be calculated. The script reports `null` for `pearson_r` in those cases.

## Preparing your own files

Supply two tables: one for the study and one for the reference panel. Each must have a header row and these columns:

| Column | Contents |
|---|---|
| `SNP` | A unique variant identifier, such as an rsID. Use consistent identifiers in both files. |
| `A1` | The allele whose frequency is reported. |
| `A2` | The other allele at that SNP. |
| `EAF` | Frequency of A1, expressed as a number between 0 and 1. |

An illustrative input table is:

```text
SNP    A1    A2    EAF
rs1    A     G     0.20
rs2    C     A     0.70
rs3    G     T     0.40
```

Separate columns with tabs or spaces. Plain-text files and files compressed with gzip (`.gz`) are accepted. Extra columns are allowed. Column names are case-sensitive; allele letters are converted to uppercase automatically.

**The frequency must refer to A1.** Minor allele frequency (MAF) tells you how common the less frequent allele is, but does not by itself identify that allele as A1. A value of `20` is also not a valid frequency: write `0.20` for 20%.

Both files must use compatible genome builds and variant identifiers. The script matches the `SNP` column; it does not convert genomic coordinates, check genome builds, or resolve different identifiers for the same variant. Obtain reference frequencies with explicit allele labels so they can be aligned correctly.

Run your comparison by replacing the example input paths with your own file paths. To save the output, append `> comparison.json` to the command; this writes the result to that file in your current directory.

If the frequency column has a different name, specify it with these options:

```text
--gwas-eaf A1Freq
--reference-eaf AF
```

These options identify the frequency columns; their values must still represent the corresponding A1 allele. The other required columns remain `SNP`, `A1` and `A2`. For example, a LDAK file using `Predictor` for its variant identifiers needs that header changed to `SNP` when preparing the input table.

## How the comparison works

1. The script finds SNP identifiers present in both files.
2. It aligns the alleles. If one file reports the frequency of A and the other reports G at an A/G SNP, it converts the second frequency using `1 − frequency`. It also handles complementary DNA-strand labels for unambiguous SNPs.
3. It excludes A/T and C/G SNPs, whose strand orientation cannot be resolved from the allele labels alone. It also excludes entries without two distinct A/C/G/T alleles and SNPs whose allele pairs cannot be matched.
4. It calculates the four summaries using the remaining SNPs, with equal weight for each SNP.

The script stops with an error for duplicate or missing SNP identifiers, missing required columns, inconsistent row lengths, or nonnumeric, nonfinite or out-of-range frequencies. It also stops if no usable SNPs overlap. Missing frequency values must be handled before running the comparison.

## What we tested with real data

The included example uses public files from the [LDAK-KVIK UK Biobank study](https://www.ldak-kvik.com/summaries/), described by [Hof and Speed (2025)](https://doi.org/10.1038/s41588-025-02286-z), and 1000 Genomes Phase 3 reference data.

The published height association file contains MAF but does not provide a labelled A1 frequency. We therefore used its accompanying prediction file, where `Centre` is the mean A1 allele count. For these autosomal SNPs, dividing `Centre` by two gives the A1 frequency. This convention is described in the [LDAK score-file documentation](https://dougspeed.com/profile-scores/). We did not independently establish that the prediction-model and association-analysis samples were identical.

We selected 200 overlapping SNPs from a small chromosome 1 interval, without selecting for frequency similarity. After excluding eight strand-ambiguous SNPs, the same 192 SNPs were compared with each reference group:

| 1000 Genomes reference group | Correlation | Mean absolute frequency difference |
|---|---:|---:|
| European (EUR) | 0.9994 | 0.0081 |
| Admixed American (AMR) | 0.9444 | 0.0655 |
| South Asian (SAS) | 0.9192 | 0.0789 |
| East Asian (EAS) | 0.7591 | 0.1381 |
| African (AFR) | 0.6279 | 0.1716 |

The European reference was closest in this example. All five results were checked against an independent calculation. Download sources, preprocessing steps and rerun commands are recorded in the example README.txt; full-precision outputs are in results.json.

## Using the results to choose a reference

Use frequency similarity alongside information about the study population and how the reference panel was constructed. When comparing several panels, prepare a common set of eligible SNPs for all comparisons. The script compares one pair of files at a time and does not enforce a common SNP set across separate runs.

Close allele frequencies provide supporting (but not definitive) evidence for reference selection.  The script supplies no ancestry assignment or automatic pass/fail threshold.

The included example demonstrates the calculation on real data from one small genomic interval. A broader assessment should use a representative set of SNPs across the genome and consider study composition and reference sample size when interpreting discrepancies.
