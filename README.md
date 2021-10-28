# **Cancer Mutations Pipeline**

Cancer driver identification based on the quantitative prediction of the effect of mutations in regulatory elements on transcription factor (TF) binding.

Input dataset (Liver Cancer-RIKEN, Japan; i.e. LIRI-JP): https://dcc.icgc.org/releases/current/Projects/LIRI-JP

### Filter mutation data and gene expression data
Code: `readtbl.R` <br />
Input: <br />
  1. Cancer data with somatic mutations from multiple donors. In this project, we use SNVs. <br />
  Using LIRI-JP project from ICGC (258 donors) with ~3.5 million SNVs, `readtbl.R` outputs the necessary elements for the pipeline:


  chromosome |  start |  end  |icgc_mutation_id | icgc_donor_id | ref | mut |
 ----------  |:------:|:-----:|:---------------:|:-------------:|:---:|:---:|
       chr1  |565238  |565238 | MU886122        | DO23508       |T    | C   |  
       chr1  |569604  |569604 | MU941619        | DO23508       |G    | A   |
       chr1  |570513  |570513 | MU952816        | DO23508       |A    | G   |
       chr1  |754535  |754535 | MU76676523      | DO23508       |G    | A   |
       chr1  |1487203 |1487203| MU958359        | DO23508       |C    | A   |
       chr1  |1487205 |1487205| MU76676613      | DO23508       |T    | C   |

  2. Expression data for the cancer mutations with the overlapping donors (from `icgc_donor_id` column) with the somatic mutations data.  

| icgc_donor_id|    gene_id | normalized_read_count  |         analysis_id|
|  ----------  |:    ------:| ----------------------:| ------------------:|  
|      DO50825 |   DYNC2LI1 |               11.654758| RK270_Cancer-rnaseq|
|      DO45221 |       CTTN |               15.685822|  RK100_Liver-rnaseq|
|      DO23538 |       RIC8A|                9.196332| RK131_Cancer-rnaseq|
|      DO23545 |LOC100505835|                0.312690| RK151_Cancer-rnaseq|
|      DO23542 |     HERPUD2|                6.738176|  RK143_Liver-rnaseq|

Also filter out mutations inside exons

## Filtering promoter data
Code: `filter_coding_promoters.R`<br />
Input: `dataset/input/raw/all_promoters_refseq.txt`

This script takes all promoters and take only promoters for the coding regions. Since promoters region were taken using RefSeq TSS as its center, the script filters for coding regions by looking for TSS with mappings to HGNC genes.

## Get trinucleotide background dist and mutation rate
Code: `trinucleotide_bg.R` <br />


gen_all_snps.R

## Getting the entropy
muteffect_entropy.R
--> want to put more weight in longer enhancer
https://pubmed.ncbi.nlm.nih.gov/21605215/

## Gene annotation
annotate_gene.R

## Predict the mutation effect
muteffect_entropy.R


## Combine prediction
combinepreds.R
analyze_combined_preds.R
combine_pe.R

## Top genes
analyze_top_genes.R
analyze_top_regions.R


# ** GSEA ANALYSIS **
gsea_analysis.R

R package: clusterProfiler
