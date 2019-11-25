LIRI-JP Simple Somatic Mutations
================

There are 258 donors with WGS data in the study. The total number of
single mutations is 3493148. Here, we use enhancers from the FANTOM
project, exons from RefSeq, and promoters inferred from RefSeq tss. The
three types of genomic regions are defined so that there is no overlap
among them. We also take 100,000 random intervals from human genome.
Each interval is 1,000 bp long.

The number of single mutations in each type of region:

| position |   count |   length | mutation\_rate |
| :------- | ------: | -------: | :------------- |
| enhancer |   10916 | 11759209 | 3.60e-06       |
| exon     |   68104 | 80790540 | 3.27e-06       |
| others   | 3377415 |       NA | NA             |
| promoter |   36694 | 46589695 | 3.05e-06       |
| random   |  109524 | 98325390 | 4.32e-06       |

For promoter regions, we played with its length:

| promoter      |  count |    length | mutation\_rate |
| :------------ | -----: | --------: | :------------- |
| l1000\_r1000  |  36694 |  46589695 | 3.05e-06       |
| l1000\_r0     |  18789 |  24139178 | 3.02e-06       |
| l500\_r500    |  16063 |  21491714 | 2.90e-06       |
| l5000\_r1000  | 125056 | 144491937 | 3.35e-06       |
| l10000\_r1000 | 230316 | 258443502 | 3.45e-06       |

The top-mutated
enhancers:

| chromosome |     start |       end | enhancer                 | count | length | mut\_rate |
| :--------- | --------: | --------: | :----------------------- | ----: | -----: | --------: |
| chr1       | 144533628 | 144534178 | chr1:144533627-144534178 |     7 |    550 | 0.0127273 |
| chr7       |  20254262 |  20255002 | chr7:20254261-20255002   |     7 |    740 | 0.0094595 |
| chr2       | 137013126 | 137013856 | chr2:137013125-137013856 |     6 |    730 | 0.0082192 |
| chr5       |  99387252 |  99388845 | chr5:99387251-99388845   |     6 |   1593 | 0.0037665 |
| chr8       |  58333216 |  58333500 | chr8:58333215-58333500   |     6 |    284 | 0.0211268 |
| chr12      |  65019059 |  65019979 | chr12:65019058-65019979  |     5 |    920 | 0.0054348 |
| chr14      |  34162041 |  34162301 | chr14:34162040-34162301  |     5 |    260 | 0.0192308 |
| chr17      |  57830595 |  57831044 | chr17:57830594-57831044  |     5 |    449 | 0.0111359 |
| chr21      |  40182921 |  40183779 | chr21:40182920-40183779  |     5 |    858 | 0.0058275 |
| chr5       |  16539356 |  16539674 | chr5:16539355-16539674   |     5 |    318 | 0.0157233 |

The tss whose enhancers harbor the greatest number of mutations:

| tss           | mut\_count | e\_count |
| :------------ | ---------: | -------: |
| NM\_003467    |         35 |       38 |
| NM\_003033    |         34 |       48 |
| NM\_173344    |         34 |       48 |
| NM\_001008540 |         33 |       35 |
| NM\_001160124 |         30 |       57 |
| NM\_001160125 |         30 |       57 |
| NM\_001300    |         30 |       57 |
| NM\_001171653 |         29 |       38 |
| NM\_014795    |         29 |       38 |
| NM\_014607    |         23 |       23 |

The genes whose enhancers harbor the greatest number of mutations:

| gene     | mut\_count | e\_count |
| :------- | ---------: | -------: |
| CXCR4    |         35 |       38 |
| ST3GAL1  |         34 |       48 |
| KLF6     |         30 |       57 |
| ZEB2     |         29 |       38 |
| PDE4DIP  |         25 |       26 |
| UBXN4    |         23 |       23 |
| FPR1     |         21 |       31 |
| IRF2BP2  |         21 |       41 |
| ITGB8    |         21 |       24 |
| SIGLEC14 |         21 |       28 |
