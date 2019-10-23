BRCA-EU Simple Somatic Mutations
================

There are 569 donors with WGS data in the study. The total number of
single mutations is 3594783.

The number of single mutations in each type of region:

| position |   count |   length | mutation\_rate |
| :------- | ------: | -------: | :------------- |
| enhancer |   13696 | 11800463 | 2.04e-06       |
| exon     |   89716 | 81028378 | 1.95e-06       |
| others   | 3440373 |       NA | NA             |
| promoter |   50998 | 46574079 | 1.92e-06       |

For promoter regions, we played with its length:

| promoter      |  count |    length | mutation\_rate |
| :------------ | -----: | --------: | :------------- |
| l1000\_r1000  |  50998 |  46574079 | 1.92e-06       |
| l1000\_r0     |  26456 |  24143431 | 1.93e-06       |
| l500\_r500    |  23960 |  21482947 | 1.96e-06       |
| l5000\_r1000  | 156641 | 144455017 | 1.91e-06       |
| l10000\_r1000 | 281102 | 258385125 | 1.91e-06       |

The top-mutated enhancers:

| enhancer                  | count | length |      rate |
| :------------------------ | ----: | -----: | --------: |
| chr1:177903071-177903284  |    12 |    213 | 0.0563380 |
| chr15:65596702-65597188   |     8 |    486 | 0.0164609 |
| chr10:30831188-30831972   |     6 |    784 | 0.0076531 |
| chr11:129512959-129513316 |     6 |    357 | 0.0168067 |
| chr12:69731953-69732266   |     6 |    313 | 0.0191693 |
| chr12:70634012-70634396   |     6 |    384 | 0.0156250 |
| chr17:38477236-38480096   |     6 |   2860 | 0.0020979 |
| chr1:160677742-160679048  |     5 |   1306 | 0.0038285 |
| chr10:30872538-30873163   |     5 |    625 | 0.0080000 |
| chr10:38202647-38203041   |     5 |    394 | 0.0126904 |

The genes whose enhancers harbor the greatest number of mutations:

| gene       | mut\_count | enh\_count |
| :--------- | ---------: | ---------: |
| SIRPA      |         46 |        120 |
| ELMO1      |         45 |         73 |
| CCRL2      |         42 |        102 |
| PFKFB3     |         40 |         95 |
| CEBPB      |         39 |         78 |
| CXCR4      |         36 |         73 |
| SNAI1      |         36 |         70 |
| NM\_199129 |         35 |         67 |
| FNDC3B     |         34 |         44 |
| BTG1       |         33 |         52 |
