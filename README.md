# **Cancer Mutations Pipeline**

Cancer driver identification based on quantitative prediction of the effect of mutations in regulatory elements on transcription factor (TF) binding. The TF binding change predictions are done using QBiC-Pred (https://github.com/vincentiusmartin/QBiC-Pred).

Input dataset:
1. any ICGC dataset. For example: Liver Cancer-RIKEN, Japan; i.e. LIRI-JP, https://dcc.icgc.org/releases/current/Projects/LIRI-JP.
2. Promoter, enhancer, exon data, e.g. from [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/). Processed RefSeq data is provided in `dataset/epdata` folder.

All input paths for the code are defined in each code's header.

### 1. Read ICGC simple somatic mutations (SSM) and expression data (EXP-S)
Code: `readtbl.R` <br />
Input: <br />
  1. ICGC simple somatic mutation (SSM) data, required fields:

|  chromosome |  start |  end  |icgc_mutation_id | icgc_donor_id | ref | mut |
| --- | --- | --- | --- | --- | --- | --- |
|       chr1  |565238  |565238 | MU886122        | DO23508       |T    | C   |  
|       chr1  |569604  |569604 | MU941619        | DO23508       |G    | A   |
|       chr1  |570513  |570513 | MU952816        | DO23508       |A    | G   |
|       chr1  |754535  |754535 | MU76676523      | DO23508       |G    | A   |
|       chr1  |1487203 |1487203| MU958359        | DO23508       |C    | A   |
|       chr1  |1487205 |1487205| MU76676613      | DO23508       |T    | C   |

  2. Expression data for the cancer mutations with the overlapping donors (from `icgc_donor_id` column) with the somatic mutations data. Required fields:

| icgc_donor_id|    gene_id | normalized_read_count  |         analysis_id|
| --- | --- | --- | --- |
|    DO50825 |   DYNC2LI1 |               11.654758| RK270_Cancer-rnaseq|
|    DO45221 |       CTTN |               15.685822|  RK100_Liver-rnaseq|
|    DO23538 |       RIC8A|                9.196332| RK131_Cancer-rnaseq|
|    DO23545 |LOC100505835|                0.312690| RK151_Cancer-rnaseq|
     DO23542 |     HERPUD2|                6.738176|  RK143_Liver-rnaseq|

3. List of all exons to filter out mutations inside exons

## Filtering promoter data
Code: `filter_coding_promoters.R`<br />
Output: `dataset/epdata/promoters_refseq_coding.txt`

This script takes all promoters and take only promoters for the coding regions. Since promoters region were taken using RefSeq TSS as its center, the script filters for coding regions by looking for TSS with mappings to HGNC genes.

## Get trinucleotide background dist and mutation rate
Code: `trinucleotide_bg.R` <br />
Output: `dataset/bg_entropy_data`

Compute the trinucleotide background mutations for the promoters or enhancers.

## Getting the entropy
Code: `muteffect_entropy.R`

This script is used for:
1. Compute the mutation effects on TF binding, set `process_type <- "qbic"`. The input used are prediction files from QBiC-Pred, available to download from http://qbic.genome.duke.edu/downloads. <br />
For this use case, the script accept command argument as an index for the file in the input directory, this is used for parallel processing using SLURM, see code `submit_job.sh`.
2. Compute Shannon entropy for the input regulatory element, set `process_type <- "entropy"`. The output files are written to `dataset/bg_entropy_data`.


## Annotate enhancer or promoter with its gene mapping
Code: `annotate_cisreg.R`

The output are regulatory element's coordinates with the genes mapped to each element: `dataset/cis_annotated`.

## Combine prediction across regulatory elements
Code: `combinepreds.R`

Some predictions data are provided, for the complete predictions please contact us.

## Gene expression analysis of the combined regulatory elements
Code: `analyze_combined_preds.R`
