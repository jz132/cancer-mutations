BRCA-FR Simple Somatic Mutations
================

There are 72 donors with WGS data in the study. The total number of
single mutations is 622660.

The number of single mutations in each type of region:

| position |  count |   length | mutation\_rate |
| :------- | -----: | -------: | :------------- |
| enhancer |   2517 | 11800463 | 2.96e-06       |
| exon     |  17289 | 81028378 | 2.96e-06       |
| others   | 592771 |       NA | NA             |
| promoter |  10083 | 46574079 | 3.01e-06       |

For promoter regions, we played with its length:

| promoter      | count |    length | mutation\_rate |
| :------------ | ----: | --------: | :------------- |
| l1000\_r1000  | 10083 |  46574079 | 3.01e-06       |
| l1000\_r0     |  5242 |  24143431 | 3.02e-06       |
| l500\_r500    |  4952 |  21482947 | 3.20e-06       |
| l5000\_r1000  | 29616 | 144455017 | 2.85e-06       |
| l10000\_r1000 | 52408 | 258385125 | 2.82e-06       |

The top-mutated enhancers:

| enhancer                  | count | length |      rate |
| :------------------------ | ----: | -----: | --------: |
| chr17:38477236-38480096   |     6 |   2860 | 0.0020979 |
| chr17:29907482-29908002   |     5 |    520 | 0.0096154 |
| chr17:45725037-45725352   |     4 |    315 | 0.0126984 |
| chr1:172852442-172852809  |     3 |    367 | 0.0081744 |
| chr10:121030438-121030814 |     3 |    376 | 0.0079787 |
| chr12:53264136-53264576   |     3 |    440 | 0.0068182 |
| chr19:47363167-47364689   |     3 |   1522 | 0.0019711 |
| chr2:204553496-204553837  |     3 |    341 | 0.0087977 |
| chr21:11184749-11185145   |     3 |    396 | 0.0075758 |
| chr3:50643029-50643593    |     3 |    564 | 0.0053191 |

The genes whose enhancers harbor the greatest number of mutations:

| gene          | mut\_count | enh\_count |
| :------------ | ---------: | ---------: |
| NM\_001145302 |         14 |         22 |
| RARA          |         14 |         18 |
| WIPF2         |         14 |         24 |
| NM\_183228    |         12 |         20 |
| TMC6          |         12 |         28 |
| CASC3         |         11 |         21 |
| ELMO1         |         11 |         73 |
| CCR7          |         10 |         19 |
| CXCR4         |         10 |         73 |
| PDE4DIP       |         10 |         49 |
