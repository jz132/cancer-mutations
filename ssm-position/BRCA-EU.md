BRCA-EU Simple Somatic Mutations
================

There are 569 donors with WGS data in the study. The total number of
single mutations is 3594783. Here, we use enhancers from the FANTOM
project, exons from RefSeq, and promoters inferred from RefSeq TSS. The
three types of genomic region are defined so that there is no overlap
among them. We also take 100,000 random intervals from human genome, the
lengths of which are 1,000 bp each. Here we compare the mutation rates
across each type of region.

## Mutation rates in different genomic regions

The number of single mutations in each type of region:

| position |   count |   length | mutation\_rate |
| :------- | ------: | -------: | :------------- |
| enhancer |   13648 | 11759211 | 2.04e-06       |
| exon     |   89381 | 80790540 | 1.94e-06       |
| others   | 3440711 |       NA | NA             |
| promoter |   51043 | 46589695 | 1.93e-06       |
| random   |  105777 | 98325390 | 1.89e-06       |

## Mutation rates for different promoter length

Promoter regions are usually defined as the parts of genome adjacent to
each TSS. The exact length to the left and right of the TSS varies from
study to study. We want to investigate if the definition of promoter
will affect the mutation rate estimation, so we’ve tried different
length and computed the corresponding rate.

| promoter      |  count |    length | mutation\_rate |
| :------------ | -----: | --------: | :------------- |
| l1000\_r1000  |  51043 |  46589695 | 1.93e-06       |
| l1000\_r0     |  32919 |  29235415 | 1.98e-06       |
| l500\_r500    |  23992 |  21491714 | 1.96e-06       |
| l5000\_r1000  | 167637 | 153418487 | 1.92e-06       |
| l10000\_r1000 | 303310 | 275723250 | 1.93e-06       |

## Exon mutations

<!-- The transcripts whose exon harbor the greatest number of mutations: -->

The genes whose exons harbor the greatest number of
mutations:

| gene    | count | exon\_length | expected\_count | var\_count | log.p\_value |
| :------ | ----: | -----------: | --------------: | ---------: | -----------: |
| PIK3CA  |   175 |         9114 |        9.842094 |   9.842050 |  148.5100134 |
| TP53    |   159 |         2602 |        1.280212 |   1.280207 |  265.9640060 |
| TTN     |   152 |       105931 |      137.328451 | 137.327885 |    0.9414998 |
| MUC4    |    65 |        27344 |       28.959522 |  28.959391 |    8.2297020 |
| MUC16   |    59 |        43900 |       28.199012 |  28.198897 |    6.5543533 |
| GRIN2B  |    47 |        30361 |        9.141684 |   9.141648 |   18.1233420 |
| SYNE1   |    46 |        28000 |       82.257174 |  82.256818 |    0.0000022 |
| RYR2    |    43 |        16470 |       33.649317 |  33.649168 |    1.1692675 |
| ZFHX4   |    43 |        13969 |        9.744523 |   9.744480 |   14.3884704 |
| MUC19   |    39 |        25005 |       92.016849 |  92.016485 |    0.0000000 |
| PKHD1L1 |    39 |        20038 |       18.704594 |  18.704508 |    4.5604280 |
| FAT3    |    38 |        19072 |        5.579035 |   5.579012 |   18.7054720 |
| MUC5B   |    38 |        17965 |       32.683719 |  32.683577 |    0.7049964 |
| SRSF1   |    38 |        10503 |        2.749531 |   2.749518 |   29.1891074 |
| CSMD1   |    37 |        14387 |       15.455492 |  15.455424 |    5.6329562 |
| KMT2C   |    37 |        16921 |       23.729191 |  23.729092 |    2.1573368 |
| OBSCN   |    37 |        25864 |       38.654957 |  38.654802 |    0.2030380 |
| PRRC2C  |    36 |        10418 |       12.484470 |  12.484421 |    7.3471680 |
| LAMA2   |    35 |         9757 |       17.416895 |  17.416827 |    3.8661570 |
| MUC17   |    35 |        14363 |        9.690609 |   9.690566 |    9.5658540 |

## Promoter mutations

<!-- The transcripts whose promoters harbor the greatest number of mutations: -->

The genes whose promoters harbor the greatest number of
mutations:

| gene    | count | promoter\_length | expected\_count | var\_count | log.p\_value |
| :------ | ----: | ---------------: | --------------: | ---------: | -----------: |
| PSMC5   |    39 |             4002 |        4.797273 |   4.797252 |    21.778942 |
| PLEKHS1 |    23 |             2001 |        1.944497 |   1.944490 |    16.577849 |
| FDPS    |    20 |             4002 |        4.459466 |   4.459447 |     7.234882 |
| FTSJ3   |    20 |             2001 |        2.386003 |   2.385992 |    11.816871 |
| H2AC12  |    19 |             2001 |        2.328891 |   2.328880 |    11.067268 |
| CD2BP2  |    18 |             4002 |        4.653940 |   4.653921 |     5.686860 |
| H2BC12  |    18 |             2001 |        2.290080 |   2.290069 |    10.268216 |
| TECR    |    18 |             4002 |        4.593619 |   4.593600 |     5.764399 |
| HINT3   |    17 |             2001 |        2.234189 |   2.234179 |     9.529216 |
| LITAF   |    17 |             4002 |        4.831085 |   4.831065 |     4.887534 |
| TAGLN2  |    17 |             4002 |        4.569284 |   4.569266 |     5.193344 |
| CD1A    |    16 |             4002 |        4.019626 |   4.019611 |     5.284342 |
| LMNB1   |    16 |             4002 |        4.964406 |   4.964384 |     4.196498 |
| UBR5    |    16 |             2001 |        2.395598 |   2.395589 |     8.225051 |
| ZNF32   |    16 |             4002 |        4.694223 |   4.694204 |     4.477145 |
| GRB7    |    15 |             4002 |        4.435194 |   4.435176 |     4.201427 |
| PDZRN4  |    15 |             4002 |        4.335672 |   4.335654 |     4.309571 |
| POTEB   |    15 |             8004 |        8.538849 |   8.538816 |     1.547157 |
| AFF3    |    14 |             4002 |        4.479360 |   4.479343 |     3.619097 |
| CMTR1   |    14 |             2001 |        2.271279 |   2.271269 |     6.868547 |

## Enhancer mutations

<!-- The transcripts whose enhancers harbor the greatest number of mutations: -->

The genes whose enhancers harbor the greatest number of
mutations:

| gene    | count | enhancer\_length | expected\_count | var\_count | log.p\_value |
| :------ | ----: | ---------------: | --------------: | ---------: | -----------: |
| CEBPB   |    39 |            27623 |        32.36307 |   32.36293 |    0.8506428 |
| SNAI1   |    36 |            24875 |        29.19075 |   29.19062 |    0.9093010 |
| TMEM189 |    35 |            24116 |        28.17080 |   28.17067 |    0.9258888 |
| BTG1    |    33 |            19056 |        20.84631 |   20.84623 |    2.0758026 |
| TNFAIP3 |    32 |            19082 |        20.83645 |   20.83637 |    1.8616050 |
| CD48    |    31 |            12444 |        13.78719 |   13.78714 |    4.3409247 |
| CD244   |    30 |            10997 |        12.22782 |   12.22777 |    4.9010358 |
| KLF6    |    30 |            22742 |        26.33709 |   26.33697 |    0.5813484 |
| SLC9A8  |    29 |            20532 |        24.04226 |   24.04215 |    0.7452319 |
| FNDC3B  |    27 |            15941 |        17.60389 |   17.60382 |    1.6517197 |
| PTPN1   |    27 |            18675 |        22.41314 |   22.41304 |    0.7183335 |
| ZMIZ1   |    26 |            19004 |        23.33803 |   23.33793 |    0.4984081 |
| IL1A    |    25 |            18405 |        20.76697 |   20.76688 |    0.6931530 |
| BCAT1   |    24 |            11023 |        12.15235 |   12.15230 |    2.7619156 |
| BCL2L14 |    24 |            19470 |        22.17009 |   22.16999 |    0.4243373 |
| CD44    |    24 |            22004 |        24.55738 |   24.55728 |    0.2427781 |
| IL1RN   |    24 |            18697 |        21.14517 |   21.14508 |    0.5301540 |
| IL36B   |    24 |            18202 |        20.56312 |   20.56303 |    0.5991416 |
| IL36G   |    24 |            18697 |        21.14517 |   21.14508 |    0.5301540 |
| PPIF    |    24 |            17475 |        21.50023 |   21.50013 |    0.4913281 |

## Regulatory mutations

<!-- The transcripts whose regulatory regions have the most significant excess of mutations: -->

The genes whose regulatory regions have the most significant excess of
mutations:

| gene     | count | promoter\_length | enhancer\_length | total\_length | expected\_count | log.p\_value |
| :------- | ----: | ---------------: | ---------------: | ------------: | --------------: | -----------: |
| PSMC5    |    39 |             4002 |               NA |          4002 |        4.797273 |    21.778942 |
| PLEKHS1  |    23 |             2001 |               NA |          2001 |        1.944497 |    16.577849 |
| FTSJ3    |    20 |             2001 |               NA |          2001 |        2.386003 |    11.816871 |
| H2AC12   |    19 |             2001 |               NA |          2001 |        2.328891 |    11.067268 |
| H2BC12   |    18 |             2001 |               NA |          2001 |        2.290080 |    10.268216 |
| HINT3    |    17 |             2001 |               NA |          2001 |        2.234189 |     9.529216 |
| H2BC7    |    14 |             2001 |               NA |          2001 |        2.107373 |     7.258196 |
| FDPS     |    20 |             4002 |               NA |          4002 |        4.459466 |     7.234882 |
| PDZRN4   |    21 |             4002 |              502 |          4504 |        4.914081 |     7.213911 |
| CMTR1    |    14 |             2001 |               NA |          2001 |        2.271279 |     6.868547 |
| MED16    |    14 |             2001 |               NA |          2001 |        2.285222 |     6.836925 |
| CCDC107  |    15 |             2001 |              397 |          2398 |        2.724898 |     6.689597 |
| UBR5     |    23 |             2001 |             3545 |          5546 |        6.345315 |     6.580562 |
| H3C8     |    13 |             2001 |               NA |          2001 |        2.264833 |     6.086801 |
| H2AC15   |    13 |             2001 |               NA |          2001 |        2.311296 |     5.990655 |
| OR4C15   |    12 |             2001 |               13 |          2014 |        1.978227 |     5.913416 |
| FAM163A  |    15 |             2001 |              589 |          2590 |        3.177076 |     5.871218 |
| TECR     |    18 |             4002 |               NA |          4002 |        4.593619 |     5.764399 |
| CD2BP2   |    18 |             4002 |               NA |          4002 |        4.653940 |     5.686860 |
| RFPL4AL1 |    12 |             2001 |               NA |          2001 |        2.090024 |     5.671167 |
