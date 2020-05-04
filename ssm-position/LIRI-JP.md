LIRI-JP Simple Somatic Mutations
================

There are 258 donors with WGS data in the study. The total number of
single mutations is 3493148. Here, we use enhancers from the FANTOM
project, exons from RefSeq, and promoters inferred from RefSeq TSS. The
three types of genomic region are defined so that there is no overlap
among them. We also take 100,000 random intervals from human genome, the
lengths of which are 1,000 bp each. Here we compare the mutation rates
across each type of region.

## Mutation rates in different genomic regions

The number of single mutations in each type of region:

| position |   count |   length | mutation\_rate |
| :------- | ------: | -------: | :------------- |
| enhancer |   10916 | 11759211 | 3.60e-06       |
| exon     |   68104 | 80790540 | 3.27e-06       |
| others   | 3377459 |       NA | NA             |
| promoter |   36669 | 46589695 | 3.05e-06       |
| random   |  109526 | 98325390 | 4.32e-06       |

## Mutation rates for different promoter length

Promoter regions are usually defined as the parts of genome adjacent to
each TSS. The exact length to the left and right of the TSS varies from
study to study. We want to investigate if the definition of promoter
will affect the mutation rate estimation, so we’ve tried different
length and computed the corresponding rate.

| promoter      |  count |    length | mutation\_rate |
| :------------ | -----: | --------: | :------------- |
| l1000\_r1000  |  36669 |  46589695 | 3.05e-06       |
| l1000\_r0     |  23187 |  29235415 | 3.07e-06       |
| l500\_r500    |  16042 |  21491714 | 2.89e-06       |
| l5000\_r1000  | 132598 | 153418487 | 3.35e-06       |
| l10000\_r1000 | 245773 | 275723250 | 3.45e-06       |

## Exon mutations

<!-- The transcripts whose exon harbor the greatest number of mutations: -->

The genes whose exons harbor the greatest number of
mutations:

| gene      | count | exon\_length | expected\_count | var\_count | log.p\_value |
| :-------- | ----: | -----------: | --------------: | ---------: | -----------: |
| TP53      |    67 |         2602 |       2.1007753 |  2.1007677 |    73.861248 |
| CTNNB1    |    74 |         4062 |       3.3307476 |  3.3307353 |    70.278272 |
| PCLO      |    57 |        17115 |      13.7720225 | 13.7719741 |    17.548953 |
| LRP1B     |    46 |        16622 |      13.7109674 | 13.7109166 |    11.241712 |
| EYS       |    34 |        10696 |       8.6614813 |  8.6614506 |    10.231440 |
| MUC19     |    52 |        25005 |      19.1415281 | 19.1414651 |     9.364756 |
| ROBO2     |    28 |         8011 |       6.6416579 |  6.6416332 |     9.233099 |
| KRTAP5-11 |    11 |         1022 |       0.7685216 |  0.7685192 |     9.164113 |
| XIRP2     |    35 |        12684 |      10.2375253 | 10.2374892 |     8.959958 |
| REG3G     |    11 |         1049 |       0.8226532 |  0.8226504 |     8.860380 |
| DCC       |    31 |        10239 |       8.3652466 |  8.3652165 |     8.821133 |
| ALB       |    15 |         2350 |       1.8708192 |  1.8708126 |     8.794940 |
| CSMD3     |    35 |        12686 |      10.4635533 | 10.4635147 |     8.722455 |
| BAGE5     |    13 |         1764 |       1.4187239 |  1.4187189 |     8.389706 |
| BAGE4     |    13 |         1765 |       1.4194054 |  1.4194004 |     8.387267 |
| AFF4      |    29 |         9601 |       7.9335997 |  7.9335707 |     8.175749 |
| BAGE2     |    13 |         1901 |       1.5274642 |  1.5274587 |     8.016274 |
| BAGE3     |    13 |         1901 |       1.5274642 |  1.5274587 |     8.016274 |
| DPP10     |    27 |         8953 |       7.4168069 |  7.4167796 |     7.630398 |
| ZNF99     |    25 |         7885 |       6.5699329 |  6.5699088 |     7.480150 |

## Promoter mutations

<!-- The transcripts whose promoters harbor the greatest number of mutations: -->

The genes whose promoters harbor the greatest number of
mutations:

| gene     | count | promoter\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| ALG1L2   |    17 |             2001 |        1.132153 |   1.132150 |    14.098260 |
| ZNF595   |    14 |             2001 |        1.172910 |   1.172907 |    10.444922 |
| ZNF718   |    14 |             2001 |        1.174072 |   1.174069 |    10.439370 |
| IFI44L   |    14 |             2001 |        1.207000 |   1.206997 |    10.284472 |
| SEMA3D   |    11 |             2001 |        1.256644 |   1.256640 |     7.007975 |
| ALB      |    11 |             2001 |        1.275272 |   1.275268 |     6.945030 |
| CHMP4C   |    10 |             2001 |        1.183195 |   1.183192 |     6.294139 |
| SLC24A5  |    10 |             2001 |        1.193497 |   1.193494 |     6.260516 |
| CDH18    |    10 |             2001 |        1.197100 |   1.197097 |     6.248833 |
| SEMA3E   |    13 |             4002 |        2.424999 |   2.424992 |     5.764790 |
| ITPRID1  |    13 |             4002 |        2.434220 |   2.434214 |     5.747030 |
| DYNAP    |    13 |             4002 |        2.516310 |   2.516303 |     5.592431 |
| TFPI2    |     9 |             2001 |        1.183428 |   1.183425 |     5.361421 |
| ITIH2    |     9 |             2001 |        1.184582 |   1.184579 |     5.358059 |
| ZPBP     |     9 |             2001 |        1.205778 |   1.205775 |     5.296925 |
| MNDA     |     9 |             2001 |        1.240076 |   1.240073 |     5.200541 |
| LINGO2   |    12 |             4002 |        2.356237 |   2.356231 |     5.151590 |
| NBPF15   |    12 |             4002 |        2.374016 |   2.374010 |     5.119434 |
| GFRAL    |     9 |             2001 |        1.290817 |   1.290813 |     5.063378 |
| SERPINC1 |     8 |             2001 |        1.124122 |   1.124119 |     4.630108 |

## Enhancer mutations

<!-- The transcripts whose enhancers harbor the greatest number of mutations: -->

The genes whose enhancers harbor the greatest number of
mutations:

| gene     | count | enhancer\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| GRIA4    |    10 |             1210 |       1.0879144 |  1.0879102 |     6.621484 |
| CXCR4    |    35 |            15721 |      14.2075841 | 14.2075279 |     5.632876 |
| PDE4DIP  |    25 |            10498 |       9.5228213 |  9.5227837 |     4.663985 |
| MMRN1    |     6 |              581 |       0.5554195 |  0.5554170 |     4.595306 |
| SATB1    |    12 |             3137 |       2.8016570 |  2.8016462 |     4.424852 |
| UBXN4    |    23 |             9620 |       8.6843974 |  8.6843632 |     4.403021 |
| MACC1    |     8 |             1355 |       1.2153779 |  1.2153732 |     4.393655 |
| RGMB     |     9 |             1769 |       1.5775059 |  1.5774999 |     4.389865 |
| TRPS1    |    15 |             4648 |       4.2748069 |  4.2747893 |     4.377381 |
| ZEB2     |    29 |            13984 |      12.5706861 | 12.5706373 |     4.295168 |
| ITGB8    |    21 |             8668 |       7.7953655 |  7.7953355 |     4.180372 |
| PXDN     |     7 |             1002 |       0.9817526 |  0.9817480 |     4.128816 |
| RBFOX1   |    19 |             7440 |       6.7192160 |  6.7191898 |     4.111110 |
| KIAA1826 |     6 |              905 |       0.8128073 |  0.8128042 |     3.697683 |
| ST3GAL1  |    34 |            18914 |      17.1774116 | 17.1773438 |     3.658291 |
| MMP16    |     5 |              656 |       0.5696099 |  0.5696078 |     3.505999 |
| ANO9     |     4 |              354 |       0.3304603 |  0.3304589 |     3.417903 |
| LRRC56   |     4 |              354 |       0.3304603 |  0.3304589 |     3.417903 |
| OR2G6    |     4 |              375 |       0.3312642 |  0.3312629 |     3.413959 |
| OR2L3    |     4 |              375 |       0.3312642 |  0.3312629 |     3.413959 |

## Regulatory mutations

<!-- The transcripts whose regulatory regions have the most significant excess of mutations: -->

The genes whose regulatory regions have the most significant excess of
mutations:

| gene    | count | promoter\_length | enhancer\_length | total\_length | expected\_count | log.p\_value |
| :------ | ----: | ---------------: | ---------------: | ------------: | --------------: | -----------: |
| ALG1L2  |    17 |             2001 |               NA |          2001 |        1.132153 |    14.098260 |
| ZNF595  |    14 |             2001 |               NA |          2001 |        1.172910 |    10.444922 |
| ZNF718  |    14 |             2001 |               NA |          2001 |        1.174072 |    10.439370 |
| IFI44L  |    14 |             2001 |               NA |          2001 |        1.207000 |    10.284472 |
| SEMA3D  |    11 |             2001 |               NA |          2001 |        1.256644 |     7.007975 |
| ALB     |    11 |             2001 |               43 |          2044 |        1.305237 |     6.845901 |
| CDH18   |    11 |             2001 |              151 |          2152 |        1.333964 |     6.753229 |
| RBFOX1  |    28 |             4002 |             7440 |         11442 |        9.073156 |     6.447168 |
| CHMP4C  |    10 |             2001 |               NA |          2001 |        1.183195 |     6.294139 |
| MNDA    |    15 |             2001 |             1979 |          3980 |        3.026905 |     6.126324 |
| MS4A4A  |    12 |             2001 |              790 |          2791 |        1.899550 |     6.093755 |
| SEMA3E  |    13 |             4002 |               NA |          4002 |        2.424999 |     5.764790 |
| ITPRID1 |    13 |             4002 |               NA |          4002 |        2.434220 |     5.747030 |
| DYNAP   |    13 |             4002 |               NA |          4002 |        2.516310 |     5.592431 |
| ZPBP    |     9 |             2001 |               NA |          2001 |        1.205778 |     5.296925 |
| TFPI2   |    10 |             2001 |              435 |          2436 |        1.574397 |     5.206245 |
| NBPF15  |    12 |             4002 |               NA |          4002 |        2.374016 |     5.119434 |
| LINGO2  |    14 |             4002 |             1082 |          5084 |        3.309372 |     4.994982 |
| CXCR4   |    35 |             2001 |            15721 |         17722 |       15.360229 |     4.925451 |
| TRIML1  |    10 |             2001 |              643 |          2644 |        1.728307 |     4.861099 |
