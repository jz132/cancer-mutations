BRCA-FR Simple Somatic Mutations
================

There are 72 donors with WGS data in the study. The total number of
single mutations is 622660. Here, we use enhancers from the FANTOM
project, exons from RefSeq, and promoters inferred from RefSeq TSS. The
three types of genomic region are defined so that there is no overlap
among them. We also take 100,000 random intervals from human genome, the
lengths of which are 1,000 bp each. Here we compare the mutation rates
across each type of region.

## Mutation rates in different genomic regions

The number of single mutations in each type of region:

| position |  count |   length | mutation\_rate |
| :------- | -----: | -------: | :------------- |
| enhancer |   2506 | 11759211 | 2.96e-06       |
| exon     |  17238 | 80790540 | 2.96e-06       |
| others   | 592836 |       NA | NA             |
| promoter |  10080 | 46589695 | 3.00e-06       |
| random   |  19598 | 98325390 | 2.77e-06       |

## Mutation rates for different promoter length

Promoter regions are usually defined as the parts of genome adjacent to
each TSS. The exact length to the left and right of the TSS varies from
study to study. We want to investigate if the definition of promoter
will affect the mutation rate estimation, so we’ve tried different
length and computed the corresponding rate.

| promoter      | count |    length | mutation\_rate |
| :------------ | ----: | --------: | :------------- |
| l1000\_r1000  | 10080 |  46589695 | 3.00e-06       |
| l1000\_r0     |  6504 |  29235415 | 3.09e-06       |
| l500\_r500    |  4951 |  21491714 | 3.20e-06       |
| l5000\_r1000  | 31720 | 153418487 | 2.87e-06       |
| l10000\_r1000 | 56498 | 275723250 | 2.85e-06       |

## Exon mutations

<!-- The transcripts whose exon harbor the greatest number of mutations: -->

The genes whose exons harbor the greatest number of
mutations:

| gene     | count | exon\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | -----------: | --------------: | ---------: | -----------: |
| TP53     |    26 |         2602 |       0.5474910 |  0.5474858 |    33.636699 |
| CDK12    |    25 |        12233 |       2.5214140 |  2.5213886 |    16.200459 |
| PIK3CA   |    21 |         9114 |       1.7084653 |  1.7084482 |    15.530615 |
| CYB5D1   |    19 |         6605 |       1.3743065 |  1.3742929 |    15.027553 |
| SPAG9    |    19 |         8273 |       1.6742530 |  1.6742363 |    13.521814 |
| SRSF1    |    18 |        10503 |       1.9846123 |  1.9845931 |    11.262462 |
| KIRREL2  |    12 |         3863 |       0.8683480 |  0.8683397 |     9.763262 |
| MED1     |    13 |         8165 |       1.6777828 |  1.6777660 |     7.546446 |
| PDZK1IP1 |     6 |          898 |       0.1974816 |  0.1974797 |     7.157555 |
| NFIA     |    13 |         9493 |       1.8293190 |  1.8293013 |     7.118796 |
| ERBB2    |    10 |         5237 |       1.1271894 |  1.1271788 |     6.482829 |
| GSDMB    |     7 |         2030 |       0.4387368 |  0.4387325 |     6.373206 |
| SP6      |     9 |         4029 |       0.8983527 |  0.8983445 |     6.328383 |
| NAA38    |     7 |         1767 |       0.4947078 |  0.4947023 |     6.029320 |
| RSBN1    |    10 |         6628 |       1.2728780 |  1.2728658 |     6.011885 |
| SLCO5A1  |    10 |         5425 |       1.3040814 |  1.3040680 |     5.918895 |
| GRB7     |     7 |         2348 |       0.5406524 |  0.5406472 |     5.776663 |
| ORMDL3   |     8 |         3799 |       0.7889204 |  0.7889132 |     5.732423 |
| XRCC3    |     7 |         2612 |       0.5720774 |  0.5720721 |     5.616752 |
| KDM6B    |    10 |         6716 |       1.5086944 |  1.5086800 |     5.365770 |

## Promoter mutations

<!-- The transcripts whose promoters harbor the greatest number of mutations: -->

The genes whose promoters harbor the greatest number of
mutations:

| gene     | count | promoter\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| TLCD1    |    10 |             4002 |       0.7266916 |  0.7266858 |     8.232355 |
| H2BC7    |     7 |             2001 |       0.3069854 |  0.3069830 |     7.409011 |
| PHB      |     8 |             4002 |       0.6593197 |  0.6593145 |     6.306314 |
| LITAF    |     8 |             4002 |       0.7013188 |  0.7013135 |     6.107843 |
| APLP1    |     6 |             2001 |       0.3152287 |  0.3152266 |     5.982587 |
| LY6G5C   |     6 |             2001 |       0.3187640 |  0.3187616 |     5.954834 |
| PECAM1   |     6 |             2001 |       0.3214091 |  0.3214065 |     5.934279 |
| UVSSA    |     8 |             4002 |       0.7476713 |  0.7476653 |     5.903223 |
| GAGE7    |     6 |             2001 |       0.3455269 |  0.3455240 |     5.754660 |
| GAGE4    |     6 |             2001 |       0.3455828 |  0.3455799 |     5.754260 |
| GAGE5    |     6 |             2001 |       0.3465466 |  0.3465437 |     5.747359 |
| PPP1R1B  |     7 |             4002 |       0.6624398 |  0.6624350 |     5.204940 |
| USP17L18 |     7 |             4002 |       0.7218232 |  0.7218168 |     4.966287 |
| PIGU     |     5 |             2001 |       0.3288652 |  0.3288626 |     4.612637 |
| TOM1L1   |     5 |             2001 |       0.3307986 |  0.3307960 |     4.600603 |
| AP2A1    |     5 |             2001 |       0.3474244 |  0.3474217 |     4.500087 |
| TSNARE1  |     5 |             2001 |       0.3494489 |  0.3494464 |     4.488196 |
| PNMT     |     5 |             2001 |       0.3522571 |  0.3522544 |     4.471823 |
| FAM102A  |     5 |             2001 |       0.3686254 |  0.3686225 |     4.379067 |
| ERN1     |     5 |             2001 |       0.3812891 |  0.3812858 |     4.310263 |

## Enhancer mutations

<!-- The transcripts whose enhancers harbor the greatest number of mutations: -->

The genes whose enhancers harbor the greatest number of
mutations:

| gene    | count | enhancer\_length | expected\_count | var\_count | log.p\_value |
| :------ | ----: | ---------------: | --------------: | ---------: | -----------: |
| SP2     |     7 |             1734 |       0.3656982 |  0.3656946 |     6.899179 |
| RARA    |    14 |            11339 |       2.3745660 |  2.3745440 |     6.639561 |
| WIPF2   |    14 |            12638 |       2.7597579 |  2.7597318 |     5.879765 |
| IKZF3   |    12 |             9498 |       2.0611848 |  2.0611652 |     5.732167 |
| ZPBP2   |     7 |             2940 |       0.6063682 |  0.6063629 |     5.452700 |
| CASC3   |    11 |            10667 |       2.3409181 |  2.3408957 |     4.462133 |
| MUSTN1  |     3 |              372 |       0.0799956 |  0.0799948 |     4.094957 |
| CCR7    |    10 |            10283 |       2.1808146 |  2.1807939 |     4.026769 |
| EVI2A   |     5 |             2153 |       0.4537302 |  0.4537258 |     3.958494 |
| COX15   |     2 |               98 |       0.0172049 |  0.0172048 |     3.834702 |
| TMC8    |     7 |             4427 |       1.1095312 |  1.1095198 |     3.804647 |
| BAGE    |     3 |              799 |       0.1241179 |  0.1241168 |     3.536950 |
| BAGE2   |     3 |              799 |       0.1241179 |  0.1241168 |     3.536950 |
| PFDN5   |     6 |             3951 |       0.9056569 |  0.9056481 |     3.449822 |
| PRKCQ   |     5 |             3159 |       0.6352387 |  0.6352324 |     3.292596 |
| TMC6    |     8 |             7727 |       1.8749873 |  1.8749680 |     3.137095 |
| GNA13   |     6 |             4805 |       1.0630065 |  1.0629955 |     3.089848 |
| ABCC2   |     2 |              237 |       0.0455996 |  0.0455991 |     2.996286 |
| GLT8D1  |     3 |              985 |       0.2058659 |  0.2058639 |     2.904104 |
| FAM107A |     3 |              941 |       0.2132756 |  0.2132736 |     2.860421 |

## Regulatory mutations

<!-- The transcripts whose regulatory regions have the most significant excess of mutations: -->

The genes whose regulatory regions have the most significant excess of
mutations:

| gene     | count | promoter\_length | enhancer\_length | total\_length | expected\_count | log.p\_value |
| :------- | ----: | ---------------: | ---------------: | ------------: | --------------: | -----------: |
| TLCD1    |    10 |             4002 |               NA |          4002 |       0.7266916 |     8.232355 |
| H2BC7    |     7 |             2001 |               NA |          2001 |       0.3069854 |     7.409011 |
| RARA     |    15 |             2001 |            11339 |         13340 |       2.6889505 |     6.761640 |
| APLP1    |     7 |             2001 |              560 |          2561 |       0.4490657 |     6.306366 |
| PHB      |     8 |             4002 |               NA |          4002 |       0.6593197 |     6.306314 |
| ZPBP2    |     9 |             2001 |             2940 |          4941 |       0.9768676 |     6.031293 |
| PPP1R1B  |     8 |             4002 |              384 |          4386 |       0.7372705 |     5.947914 |
| UVSSA    |     8 |             4002 |               NA |          4002 |       0.7476713 |     5.903223 |
| PECAM1   |     9 |             2001 |             3386 |          5387 |       1.0212055 |     5.874958 |
| IKZF3    |    13 |             2001 |             9498 |         11499 |       2.3897185 |     5.833490 |
| GAGE7    |     6 |             2001 |               NA |          2001 |       0.3455269 |     5.754660 |
| GAGE4    |     6 |             2001 |               NA |          2001 |       0.3455828 |     5.754260 |
| GAGE5    |     6 |             2001 |               NA |          2001 |       0.3465466 |     5.747359 |
| TOM1L1   |     6 |             2001 |              226 |          2227 |       0.3716406 |     5.574469 |
| WIPF2    |    14 |             2001 |            12638 |         14639 |       3.0953425 |     5.316120 |
| PFDN5    |     9 |             2001 |             3951 |          5952 |       1.2418225 |     5.195714 |
| USP17L18 |     7 |             4002 |               NA |          4002 |       0.7218232 |     4.966287 |
| SP2      |     7 |             2001 |             1734 |          3735 |       0.7250411 |     4.953974 |
| LY6G5C   |     6 |             2001 |              836 |          2837 |       0.4974767 |     4.861039 |
| ACSL4    |     8 |             4002 |             1774 |          5776 |       1.1038113 |     4.685721 |
