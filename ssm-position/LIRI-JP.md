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

| gene    | count | exon\_length | expected\_count |  var\_count | log.p\_value |
| :------ | ----: | -----------: | --------------: | ----------: | -----------: |
| TTN     |   113 |       105931 |     101.6879587 | 101.6875842 |    0.8469392 |
| CTNNB1  |    74 |         4062 |       3.6663912 |   3.6663775 |   67.3364204 |
| TP53    |    67 |         2602 |       0.9525031 |   0.9524997 |   96.3854415 |
| PCLO    |    57 |        17115 |       5.4078092 |   5.4077894 |   37.1316781 |
| MUC16   |    52 |        43900 |      21.5072022 |  21.5071235 |    7.7301925 |
| MUC19   |    52 |        25005 |      70.7616736 |  70.7614189 |    0.0037072 |
| LRP1B   |    46 |        16622 |      21.4504388 |  21.4503557 |    5.5514074 |
| CSMD3   |    35 |        12686 |      13.4397688 |  13.4397152 |    6.1582298 |
| XIRP2   |    35 |        12684 |       3.0543803 |   3.0543682 |   24.3300191 |
| EYS     |    34 |        10696 |      12.8349588 |  12.8349092 |    6.1642473 |
| CSMD1   |    33 |        14387 |      10.7295236 |  10.7294811 |    7.4271114 |
| ADGRV1  |    31 |        19423 |      40.0911829 |  40.0910357 |    0.0268858 |
| CCDC141 |    31 |        14275 |       4.0035270 |   4.0035123 |   16.9202235 |
| DCC     |    31 |        10239 |       7.6432366 |   7.6432025 |    9.7355531 |
| GUCY1A2 |    31 |        16228 |       2.6913362 |   2.6913249 |   21.7168034 |
| NBEA    |    31 |        13106 |      24.0532501 |  24.0531598 |    1.0096926 |
| ABCA13  |    30 |        17246 |      17.0809349 |  17.0808697 |    2.5353476 |
| ADGRL3  |    30 |        12664 |       5.8631880 |   5.8631666 |   11.8356598 |
| DST     |    30 |        24034 |      30.7393969 |  30.7392770 |    0.2386889 |
| PIAS2   |    30 |        20838 |       7.3125397 |   7.3125102 |    9.5618051 |

## Promoter mutations

<!-- The transcripts whose promoters harbor the greatest number of mutations: -->

The genes whose promoters harbor the greatest number of
mutations:

| gene    | count | promoter\_length | expected\_count | var\_count | log.p\_value |
| :------ | ----: | ---------------: | --------------: | ---------: | -----------: |
| ALG1L2  |    17 |             2001 |        1.419025 |   1.419021 |    12.548016 |
| IFI44L  |    17 |             2001 |        1.477468 |   1.477463 |    12.273901 |
| LINGO2  |    16 |             4002 |        2.939189 |   2.939180 |     7.024020 |
| ZNF595  |    16 |             2001 |        1.480872 |   1.480867 |    11.196108 |
| ZNF718  |    16 |             2001 |        1.482730 |   1.482725 |    11.188151 |
| CD24    |    15 |             4002 |        2.984518 |   2.984509 |     6.201158 |
| SEMA3E  |    15 |             4002 |        2.967945 |   2.967936 |     6.230771 |
| FGG     |    14 |             4002 |        2.977660 |   2.977650 |     5.504800 |
| ITPRID1 |    14 |             4002 |        2.985346 |   2.985337 |     5.492196 |
| NBPF15  |    14 |             4002 |        2.930378 |   2.930369 |     5.583231 |
| OR5AS1  |    14 |             2001 |        1.531695 |   1.531690 |     8.966765 |
| ALB     |    13 |             2001 |        1.542048 |   1.542043 |     7.968459 |
| DYNAP   |    13 |             4002 |        3.053915 |   3.053905 |     4.712659 |
| OGFOD2  |    13 |             4002 |        2.938328 |   2.938319 |     4.884671 |
| DMD     |    12 |             4002 |        3.062029 |   3.062019 |     4.064097 |
| FCRL3   |    12 |             4002 |        2.961264 |   2.961255 |     4.198887 |
| GCSAML  |    12 |             4002 |        2.986817 |   2.986808 |     4.164154 |
| SEMA3D  |    12 |             2001 |        1.521479 |   1.521474 |     7.100384 |
| AKR1C2  |    11 |             4002 |        2.966967 |   2.966958 |     3.574163 |
| CAPSL   |    11 |             4002 |        2.930485 |   2.930476 |     3.619075 |

## Enhancer mutations

<!-- The transcripts whose enhancers harbor the greatest number of mutations: -->

The genes whose enhancers harbor the greatest number of
mutations:

| gene     | count | enhancer\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| CXCR4    |    35 |            15721 |       14.207584 |  14.207528 |    5.6328757 |
| ST3GAL1  |    34 |            18914 |       17.177412 |  17.177344 |    3.6582912 |
| KLF6     |    30 |            22742 |       20.860135 |  20.860050 |    1.4579454 |
| ZEB2     |    29 |            13984 |       12.570686 |  12.570637 |    4.2951676 |
| PDE4DIP  |    25 |            10498 |        9.522821 |   9.522784 |    4.6639846 |
| UBXN4    |    23 |             9620 |        8.684397 |   8.684363 |    4.4030211 |
| FPR1     |    21 |            11798 |       10.867978 |  10.867933 |    2.3879164 |
| IRF2BP2  |    21 |            15256 |       14.046949 |  14.046891 |    1.3074655 |
| ITGB8    |    21 |             8668 |        7.795365 |   7.795336 |    4.1803723 |
| SIGLEC14 |    21 |            10942 |       10.075980 |  10.075939 |    2.7605646 |
| CD44     |    20 |            22004 |       19.781696 |  19.781620 |    0.2922163 |
| CYTH1    |    20 |            23018 |       22.147104 |  22.147005 |    0.1520080 |
| ETS1     |    20 |            10829 |        9.877970 |   9.877930 |    2.5195298 |
| FLI1     |    20 |            12184 |       11.122011 |  11.121966 |    1.9846236 |
| CEBPB    |    19 |            27623 |       25.457921 |  25.457817 |    0.0354988 |
| NR3C1    |    19 |            17368 |       15.697665 |  15.697603 |    0.6325509 |
| RBFOX1   |    19 |             7440 |        6.719216 |   6.719190 |    4.1111097 |
| BTG1     |    18 |            19056 |       17.114366 |  17.114300 |    0.3497048 |
| CCR1     |    18 |            18394 |       16.446294 |  16.446231 |    0.4169984 |
| FPR2     |    18 |            10836 |       10.006438 |  10.006397 |    1.8428462 |

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
