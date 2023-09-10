# Bam2Bsseq

The goal of Bam2Bsseq is to read in nanopore data into Bsseq objects.

## Part 1: Convert modified BAM files to BED files

Use modbam2bed program developed by Oxford Nanopore Technologies to convert modified-base BAM files to bedMethyl files.
See <https://github.com/epi2me-labs/modbam2bed> for details.

## Part 2: Convert BED files into bsseq objects

Use `Bed2Bsseq()` function to read in, extract methylation, coverage and ambiguous modification status data and create Bsseq object from bedMethyl files obtained by modbam2bed program. The resulting Bsseq object can be used for further differential methylation analysis.

## Installation

You can install the development version of Bam2Bsseq from [GitHub](https://github.com/yixuan-chen-elisa/Bam2Bsseq) with:

``` r
# install.packages("Bam2Bsseq")
devtools::install_github(repo = "yixuan-chen-elisa/Bam2Bsseq")
```
