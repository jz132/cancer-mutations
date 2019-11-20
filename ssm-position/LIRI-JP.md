LIRI-JP Simple Somatic Mutations
================

There are 258 donors with WGS data in the study. The total number of
single mutations is 3493148.

The number of single mutations in each type of region:

| position |   count |   length | mutation\_rate |
| :------- | ------: | -------: | :------------- |
| enhancer |   10953 | 11800463 | 3.60e-06       |
| exon     |   68300 | 81028378 | 3.27e-06       |
| others   | 3377186 |       NA | NA             |
| promoter |   36690 | 46574079 | 3.05e-06       |
| random   |  109626 | 98422278 | 4.32e-06       |

For promoter regions, we played with its length:

| promoter      |  count |    length | mutation\_rate |
| :------------ | -----: | --------: | :------------- |
| l1000\_r1000  |  36690 |  46574079 | 3.05e-06       |
| l1000\_r0     |  18787 |  24143431 | 3.02e-06       |
| l500\_r500    |  16060 |  21482947 | 2.90e-06       |
| l5000\_r1000  | 125042 | 144455017 | 3.36e-06       |
| l10000\_r1000 | 230276 | 258385125 | 3.45e-06       |

The top-mutated enhancers:

| enhancer                 | count | length | mut\_rate |
| :----------------------- | ----: | -----: | --------: |
| chr1:144533627-144534178 |     7 |    551 | 0.0127042 |
| chr7:20254261-20255002   |     7 |    741 | 0.0094467 |
| chr2:137013125-137013856 |     6 |    731 | 0.0082079 |
| chr5:99387251-99388845   |     6 |   1594 | 0.0037641 |
| chr8:58333215-58333500   |     6 |    285 | 0.0210526 |
| chr12:65019058-65019979  |     5 |    921 | 0.0054289 |
| chr14:34162040-34162301  |     5 |    261 | 0.0191571 |
| chr17:57830594-57831044  |     5 |    450 | 0.0111111 |
| chr21:40182920-40183779  |     5 |    859 | 0.0058207 |
| chr5:16539355-16539674   |     5 |    319 | 0.0156740 |

The genes whose enhancers harbor the greatest number of mutations:

| gene          | mut\_count | enh\_count |
| :------------ | ---------: | ---------: |
| CXCR4         |         68 |         73 |
| FLI1          |         40 |         70 |
| PDE4DIP       |         38 |         49 |
| RNF145        |         38 |         60 |
| CCRL2         |         34 |        102 |
| NM\_173344    |         34 |         48 |
| BRD2          |         30 |         51 |
| NM\_001160125 |         30 |         57 |
| PFKFB3        |         30 |         95 |
| NM\_014795    |         29 |         38 |
