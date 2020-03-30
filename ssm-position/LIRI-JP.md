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
| others   | 3377433 |       NA | NA             |
| promoter |   36695 | 46589695 | 3.05e-06       |
| random   |  109526 | 98325390 | 4.32e-06       |

## Mutation rates for different promoter length

Promoter regions are usually defined as the parts of genome adjacent to
each TSS. The exact length to the left and right of the TSS varies from
study to study. We want to investigate if the definition of promoter
will affect the mutation rate estimation, so we’ve tried different
length and computed the corresponding rate.

| promoter      |  count |    length | mutation\_rate |
| :------------ | -----: | --------: | :------------- |
| l1000\_r1000  |  36695 |  46589695 | 3.05e-06       |
| l1000\_r0     |  23198 |  29235415 | 3.08e-06       |
| l500\_r500    |  16064 |  21491714 | 2.90e-06       |
| l5000\_r1000  | 132632 | 153418487 | 3.35e-06       |
| l10000\_r1000 | 245818 | 275723250 | 3.46e-06       |

## Exon mutations

<!-- The transcripts whose exon harbor the greatest number of mutations: -->

The genes whose exons harbor the greatest number of
mutations:

| gene    | count | exon\_length | expected\_count | var\_count | log.p\_value |
| :------ | ----: | -----------: | --------------: | ---------: | -----------: |
| TP53    |    67 |         2602 |       0.9525031 |  0.9524997 |     96.38544 |
| CTNNB1  |    74 |         4062 |       3.6663912 |  3.6663775 |     67.33642 |
| SOX11   |    26 |         8720 |       0.1737769 |  0.1737762 |     46.43849 |
| NPAP1   |    23 |         8054 |       0.0964076 |  0.0964072 |     45.81806 |
| PCLO    |    57 |        17115 |       5.4078092 |  5.4077894 |     37.13168 |
| SLITRK1 |    25 |         9938 |       0.3652201 |  0.3652188 |     36.27925 |
| ZNF99   |    25 |         7885 |       0.4771562 |  0.4771542 |     33.42332 |
| PYGO1   |    27 |        16034 |       0.6500282 |  0.6500257 |     33.35992 |
| LSM8    |    19 |        12524 |       0.2413283 |  0.2413274 |     28.91508 |
| ZNF727  |    20 |         8494 |       0.3753680 |  0.3753665 |     27.05217 |
| ZNF208  |    21 |        10832 |       0.4843598 |  0.4843580 |     26.52051 |
| PCDH17  |    24 |         8013 |       0.8116004 |  0.8115970 |     26.30665 |
| PCDH10  |    24 |        13016 |       0.8131324 |  0.8131294 |     26.28763 |
| TRPS1   |    25 |        10164 |       1.0077527 |  1.0077488 |     25.52732 |
| CCDC168 |    16 |        21474 |       0.1784057 |  0.1784050 |     25.37098 |
| SLC5A3  |    15 |        11625 |       0.1344962 |  0.1344957 |     25.24059 |
| XIRP2   |    35 |        12684 |       3.0543803 |  3.0543682 |     24.33002 |
| CD24    |    14 |         4714 |       0.1118408 |  0.1118403 |     24.30533 |
| FLG2    |    17 |         9125 |       0.2780967 |  0.2780956 |     24.11376 |
| GPR26   |    17 |        10307 |       0.2873124 |  0.2873113 |     23.87684 |

## Promoter mutations

<!-- The transcripts whose promoters harbor the greatest number of mutations: -->

The genes whose promoters harbor the greatest number of
mutations:

| gene      | count | promoter\_length | expected\_count | var\_count | log.p\_value |
| :-------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| ALG1L2    |    17 |             2001 |        1.419025 |   1.419021 |    12.548016 |
| IFI44L    |    17 |             2001 |        1.477468 |   1.477463 |    12.273901 |
| ZNF595    |    16 |             2001 |        1.480872 |   1.480867 |    11.196108 |
| ZNF718    |    16 |             2001 |        1.482730 |   1.482725 |    11.188151 |
| OR5AS1    |    14 |             2001 |        1.531695 |   1.531690 |     8.966765 |
| ALB       |    13 |             2001 |        1.542048 |   1.542043 |     7.968459 |
| SEMA3D    |    12 |             2001 |        1.521479 |   1.521474 |     7.100384 |
| LINGO2    |    16 |             4002 |        2.939189 |   2.939180 |     7.024020 |
| KRTAP5-11 |    11 |             2001 |        1.403841 |   1.403836 |     6.536871 |
| SLC24A5   |    11 |             2001 |        1.473242 |   1.473237 |     6.333703 |
| CDH18     |    11 |             2001 |        1.485264 |   1.485259 |     6.299614 |
| GABRG2    |    11 |             2001 |        1.489454 |   1.489449 |     6.287807 |
| H2AC11    |    11 |             2001 |        1.495127 |   1.495123 |     6.271879 |
| SEMA3E    |    15 |             4002 |        2.967945 |   2.967936 |     6.230771 |
| CD24      |    15 |             4002 |        2.984518 |   2.984509 |     6.201158 |
| OR4K5     |    11 |             2001 |        1.520939 |   1.520934 |     6.200277 |
| OR8H3     |    11 |             2001 |        1.554463 |   1.554458 |     6.109324 |
| NBPF15    |    14 |             4002 |        2.930378 |   2.930369 |     5.583231 |
| FGG       |    14 |             4002 |        2.977660 |   2.977650 |     5.504800 |
| ITPRID1   |    14 |             4002 |        2.985346 |   2.985337 |     5.492196 |

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

| gene      | count | promoter\_length | enhancer\_length | total\_length | expected\_count | log.p\_value |
| :-------- | ----: | ---------------: | ---------------: | ------------: | --------------: | -----------: |
| ALG1L2    |    17 |             2001 |               NA |          2001 |        1.419025 |    12.548016 |
| IFI44L    |    17 |             2001 |               NA |          2001 |        1.477468 |    12.273901 |
| ZNF595    |    16 |             2001 |               NA |          2001 |        1.480872 |    11.196108 |
| ZNF718    |    16 |             2001 |               NA |          2001 |        1.482730 |    11.188151 |
| OR5AS1    |    14 |             2001 |               NA |          2001 |        1.531695 |     8.966765 |
| ALB       |    13 |             2001 |               43 |          2044 |        1.572013 |     7.871789 |
| SEMA3D    |    12 |             2001 |               NA |          2001 |        1.521479 |     7.100384 |
| OR2L3     |    13 |             2001 |              375 |          2376 |        1.847915 |     7.069120 |
| RBFOX1    |    30 |             4002 |             7440 |         11442 |        9.677225 |     6.894001 |
| CDH18     |    12 |             2001 |              151 |          2152 |        1.622128 |     6.806521 |
| LINGO2    |    18 |             4002 |             1082 |          5084 |        3.892324 |     6.774804 |
| OR4K5     |    12 |             2001 |              237 |          2238 |        1.717066 |     6.547772 |
| OR2T33    |    12 |             2001 |              375 |          2376 |        1.822949 |     6.277911 |
| H2AC11    |    11 |             2001 |               NA |          2001 |        1.495127 |     6.271879 |
| OR2M7     |    12 |             2001 |              375 |          2376 |        1.833429 |     6.252190 |
| SEMA3E    |    15 |             4002 |               NA |          4002 |        2.967945 |     6.230771 |
| KRTAP5-11 |    11 |             2001 |              126 |          2127 |        1.513201 |     6.221596 |
| CD24      |    15 |             4002 |               NA |          4002 |        2.984518 |     6.201158 |
| OR8H3     |    11 |             2001 |               NA |          2001 |        1.554463 |     6.109324 |
| GABRG2    |    11 |             2001 |              187 |          2188 |        1.657345 |     5.843658 |
