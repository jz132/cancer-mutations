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
| TP53    |   159 |         2602 |       1.2802122 |  1.2802068 |    265.96401 |
| PIK3CA  |   175 |         9114 |       9.8420935 |  9.8420500 |    148.51001 |
| CCDC168 |    31 |        21474 |       0.2466892 |  0.2466881 |     52.86214 |
| FOXA1   |    27 |         3382 |       0.3165997 |  0.3165985 |     41.65576 |
| CALN1   |    29 |         9989 |       0.6701275 |  0.6701242 |     36.26920 |
| H2AC12  |    15 |          482 |       0.0293797 |  0.0293796 |     35.10775 |
| GREM2   |    19 |         4185 |       0.1782001 |  0.1781994 |     31.39135 |
| ZNF555  |    28 |        16644 |       0.9593910 |  0.9593869 |     30.39033 |
| KCNV1   |    22 |         6983 |       0.4195467 |  0.4195450 |     29.52382 |
| CYB5D1  |    25 |         6605 |       0.6894109 |  0.6894076 |     29.51644 |
| SRSF1   |    38 |        10503 |       2.7495315 |  2.7495183 |     29.18911 |
| ZNF208  |    24 |        10832 |       0.6531753 |  0.6531727 |     28.50418 |
| ZBTB20  |    31 |        27338 |       1.6630322 |  1.6630253 |     27.76621 |
| HDGFL1  |    14 |         2215 |       0.0716536 |  0.0716533 |     26.99612 |
| TRPS1   |    29 |        10164 |       1.4398426 |  1.4398363 |     26.95939 |
| GABRA4  |    26 |        11982 |       1.0594466 |  1.0594422 |     26.39632 |
| NECTIN3 |    24 |         9161 |       0.8612724 |  0.8612682 |     25.70817 |
| TCHH    |    23 |         6998 |       0.7516071 |  0.7516030 |     25.57732 |
| FOXO3   |    19 |         7344 |       0.3840403 |  0.3840388 |     25.14031 |
| VGLL3   |    27 |        11583 |       1.3723821 |  1.3723774 |     24.89939 |

## Promoter mutations

<!-- The transcripts whose promoters harbor the greatest number of mutations: -->

The genes whose promoters harbor the greatest number of
mutations:

| gene    | count | promoter\_length | expected\_count | var\_count | log.p\_value |
| :------ | ----: | ---------------: | --------------: | ---------: | -----------: |
| PSMC5   |    39 |             4002 |        4.797273 |   4.797252 |    21.778942 |
| PLEKHS1 |    23 |             2001 |        1.944497 |   1.944490 |    16.577849 |
| FTSJ3   |    20 |             2001 |        2.386003 |   2.385992 |    11.816871 |
| H2AC12  |    19 |             2001 |        2.328891 |   2.328880 |    11.067268 |
| H2BC12  |    18 |             2001 |        2.290080 |   2.290069 |    10.268216 |
| HINT3   |    17 |             2001 |        2.234189 |   2.234179 |     9.529216 |
| UBR5    |    16 |             2001 |        2.395598 |   2.395589 |     8.225051 |
| H2BC7   |    14 |             2001 |        2.107373 |   2.107364 |     7.258196 |
| FDPS    |    20 |             4002 |        4.459466 |   4.459447 |     7.234882 |
| CMTR1   |    14 |             2001 |        2.271279 |   2.271269 |     6.868547 |
| MED16   |    14 |             2001 |        2.285222 |   2.285213 |     6.836925 |
| TOP2A   |    14 |             2001 |        2.308752 |   2.308742 |     6.784076 |
| H3C8    |    13 |             2001 |        2.264833 |   2.264823 |     6.086801 |
| CCDC107 |    13 |             2001 |        2.286994 |   2.286985 |     6.040652 |
| H2AC15  |    13 |             2001 |        2.311296 |   2.311285 |     5.990655 |
| ZNF143  |    13 |             2001 |        2.327723 |   2.327713 |     5.957211 |
| OR4C15  |    12 |             2001 |        1.972704 |   1.972696 |     5.925799 |
| TECR    |    18 |             4002 |        4.593619 |   4.593600 |     5.764399 |
| DUSP22  |    13 |             2001 |        2.437278 |   2.437268 |     5.741159 |
| CD2BP2  |    18 |             4002 |        4.653940 |   4.653921 |     5.686860 |

## Enhancer mutations

<!-- The transcripts whose enhancers harbor the greatest number of mutations: -->

The genes whose enhancers harbor the greatest number of
mutations:

| gene     | count | enhancer\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| ABCC9    |     3 |               38 |       0.0217845 |  0.0217844 |     5.770803 |
| CD244    |    30 |            10997 |      12.2278161 | 12.2277668 |     4.901036 |
| ATP10A   |     9 |             1249 |       1.5071976 |  1.5071905 |     4.541019 |
| PDZRN4   |     6 |              502 |       0.5784082 |  0.5784056 |     4.498093 |
| CD48     |    31 |            12444 |      13.7871928 | 13.7871371 |     4.340925 |
| KCNMB4   |     7 |              898 |       0.9242367 |  0.9242333 |     4.290803 |
| ADCY8    |     5 |              390 |       0.4351974 |  0.4351956 |     4.042417 |
| SHISA9   |     6 |              717 |       0.7701079 |  0.7701051 |     3.822644 |
| MYC      |    16 |             4959 |       5.4182616 |  5.4182402 |     3.770027 |
| DENND4A  |     9 |             1712 |       1.9525261 |  1.9525172 |     3.700157 |
| RAB3IP   |     9 |             1897 |       2.0582195 |  2.0582108 |     3.534550 |
| UBE3A    |    15 |             5077 |       5.1928117 |  5.1927921 |     3.475489 |
| BARX2    |     9 |             1854 |       2.2111459 |  2.2111354 |     3.312848 |
| ARHGAP30 |    21 |             7674 |       9.1959184 |  9.1958765 |     3.239923 |
| FCER1G   |    23 |             8884 |      10.5792770 | 10.5792297 |     3.201368 |
| IGHMBP2  |     9 |             1890 |       2.3556502 |  2.3556389 |     3.120515 |
| MUSTN1   |     4 |              372 |       0.4124988 |  0.4124971 |     3.060829 |
| PAPPA    |    10 |             2452 |       2.9192651 |  2.9192518 |     3.045098 |
| PRKCQ    |    11 |             3159 |       3.4596319 |  3.4596177 |     3.031454 |
| ROD1     |     5 |              633 |       0.7467208 |  0.7467173 |     2.981133 |

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
