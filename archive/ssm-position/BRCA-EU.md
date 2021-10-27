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
| others   | 3440774 |       NA | NA             |
| promoter |   50980 | 46589695 | 1.92e-06       |
| random   |  105777 | 98325390 | 1.89e-06       |

## Mutation rates for different promoter length

Promoter regions are usually defined as the parts of genome adjacent to
each TSS. The exact length to the left and right of the TSS varies from
study to study. We want to investigate if the definition of promoter
will affect the mutation rate estimation, so we’ve tried different
length and computed the corresponding rate.

| promoter      |  count |    length | mutation\_rate |
| :------------ | -----: | --------: | :------------- |
| l1000\_r1000  |  50980 |  46589695 | 1.92e-06       |
| l1000\_r0     |  32868 |  29235415 | 1.98e-06       |
| l500\_r500    |  23935 |  21491714 | 1.96e-06       |
| l5000\_r1000  | 167555 | 153418487 | 1.92e-06       |
| l10000\_r1000 | 303213 | 275723250 | 1.93e-06       |

## Exon mutations

<!-- The transcripts whose exon harbor the greatest number of mutations: -->

The genes whose exons harbor the greatest number of
mutations:

| gene   | count | exon\_length | expected\_count | var\_count | log.p\_value |
| :----- | ----: | -----------: | --------------: | ---------: | -----------: |
| TP53   |   159 |         2602 |       2.8565371 |  2.8565254 |   211.223554 |
| PIK3CA |   175 |         9114 |       8.7925087 |  8.7924751 |   156.627486 |
| H2AC12 |    15 |          482 |       0.6318011 |  0.6317981 |    15.364727 |
| FOXA1  |    27 |         3382 |       3.9807463 |  3.9807288 |    13.500561 |
| H3C4   |    15 |          923 |       1.2893227 |  1.2893160 |    10.984726 |
| SRSF1  |    38 |        10503 |      10.4292165 | 10.4291761 |    10.420553 |
| PRRC2C |    36 |        10418 |      10.9498867 | 10.9498421 |     8.756790 |
| ZFHX4  |    43 |        13969 |      14.8675832 | 14.8675215 |     8.655023 |
| LAMA2  |    35 |         9757 |      10.6583734 | 10.6583274 |     8.523425 |
| AKT1   |    19 |         3263 |       4.1995811 |  4.1995616 |     6.967060 |
| CYB5D1 |    25 |         6605 |       7.0378549 |  7.0378267 |     6.926139 |
| CALM2  |    11 |         1377 |       1.3738717 |  1.3738664 |     6.628144 |
| SPTA1  |    28 |         8068 |       8.9410215 |  8.9409828 |     6.570936 |
| GATA3  |    17 |         3076 |       3.7558720 |  3.7558546 |     6.312003 |
| GPR149 |    20 |         5055 |       5.1414912 |  5.1414710 |     6.277255 |
| H3C8   |     8 |          468 |       0.6740508 |  0.6740474 |     6.235183 |
| GREM2  |    19 |         4185 |       4.7741378 |  4.7741173 |     6.142933 |
| RYR2   |    43 |        16470 |      18.4562440 | 18.4561625 |     6.120900 |
| MUC4   |    65 |        27344 |      33.4330817 | 33.4329303 |     6.064150 |
| TRPS1  |    29 |        10164 |      10.4039332 | 10.4038922 |     5.784872 |

## Promoter mutations

<!-- The transcripts whose promoters harbor the greatest number of mutations: -->

The genes whose promoters harbor the greatest number of
mutations:

| gene    | count | promoter\_length | expected\_count | var\_count | log.p\_value |
| :------ | ----: | ---------------: | --------------: | ---------: | -----------: |
| PLEKHS1 |    22 |             2001 |        1.512917 |   1.512913 |    17.722423 |
| PSMC5   |    26 |             4002 |        3.519430 |   3.519419 |    13.865492 |
| HINT3   |    16 |             2001 |        1.663740 |   1.663735 |    10.461354 |
| UBR5    |    15 |             2001 |        1.725212 |   1.725207 |     9.263926 |
| CMTR1   |    14 |             2001 |        1.695855 |   1.695850 |     8.413790 |
| FTSJ3   |    14 |             2001 |        1.756880 |   1.756874 |     8.223385 |
| TOP2A   |    13 |             2001 |        1.700288 |   1.700283 |     7.480214 |
| DUSP22  |    13 |             2001 |        1.760178 |   1.760173 |     7.308703 |
| CD2BP2  |    17 |             4002 |        3.415922 |   3.415911 |     6.874749 |
| APLP1   |    12 |             2001 |        1.606393 |   1.606388 |     6.851075 |
| LITAF   |    17 |             4002 |        3.466268 |   3.466258 |     6.787132 |
| MRS2    |    12 |             2001 |        1.641306 |   1.641301 |     6.752879 |
| POLR2J2 |    12 |             2001 |        1.695633 |   1.695628 |     6.604732 |
| PRRC2A  |    12 |             2001 |        1.709425 |   1.709420 |     6.567985 |
| TAGLN2  |    16 |             4002 |        3.355396 |   3.355386 |     6.271841 |
| TECR    |    16 |             4002 |        3.371970 |   3.371960 |     6.244289 |
| H2BC7   |    11 |             2001 |        1.577841 |   1.577837 |     6.047216 |
| ZNF705G |    11 |             2001 |        1.586515 |   1.586510 |     6.024440 |
| MED16   |    11 |             2001 |        1.661172 |   1.661167 |     5.834145 |
| GSTA4   |    11 |             2001 |        1.696181 |   1.696176 |     5.748282 |

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

| gene    | count | promoter\_length | enhancer\_length | total\_length | expected\_count | log.p\_value |
| :------ | ----: | ---------------: | ---------------: | ------------: | --------------: | -----------: |
| PLEKHS1 |    22 |             2001 |               NA |          2001 |        1.512917 |    17.722423 |
| PSMC5   |    26 |             4002 |               NA |          4002 |        3.519430 |    13.865492 |
| HINT3   |    16 |             2001 |               NA |          2001 |        1.663740 |    10.461354 |
| CMTR1   |    14 |             2001 |               NA |          2001 |        1.695855 |     8.413790 |
| FTSJ3   |    14 |             2001 |               NA |          2001 |        1.756880 |     8.223385 |
| CD2BP2  |    17 |             4002 |               NA |          4002 |        3.415922 |     6.874749 |
| UBR5    |    22 |             2001 |             3545 |          5546 |        5.674929 |     6.807005 |
| MRS2    |    12 |             2001 |               NA |          2001 |        1.641306 |     6.752879 |
| POLR2J2 |    12 |             2001 |               NA |          2001 |        1.695633 |     6.604732 |
| APLP1   |    14 |             2001 |              560 |          2561 |        2.393694 |     6.598446 |
| PRRC2A  |    12 |             2001 |               NA |          2001 |        1.709425 |     6.567985 |
| TECR    |    16 |             4002 |               NA |          4002 |        3.371970 |     6.244289 |
| H2BC7   |    11 |             2001 |               NA |          2001 |        1.577841 |     6.047216 |
| ZNF705G |    11 |             2001 |               NA |          2001 |        1.586515 |     6.024440 |
| MED16   |    11 |             2001 |               NA |          2001 |        1.661172 |     5.834145 |
| GSTA4   |    11 |             2001 |               NA |          2001 |        1.696181 |     5.748282 |
| TPTE    |    11 |             2001 |               NA |          2001 |        1.739877 |     5.643954 |
| CCDC107 |    12 |             2001 |              397 |          2398 |        2.112849 |     5.623592 |
| SLC35E1 |    11 |             2001 |               NA |          2001 |        1.784653 |     5.540169 |
| PDZRN4  |    16 |             4002 |              502 |          4504 |        3.844206 |     5.523825 |
