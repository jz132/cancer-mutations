SKCA-BR Simple Somatic Mutations
================

There are 100 donors with WGS data in the study. The total number of
single mutations is 7593091. Here, we use enhancers from the FANTOM
project, exons from RefSeq, and promoters inferred from RefSeq TSS. The
three types of genomic region are defined so that there is no overlap
among them. We also take 100,000 random intervals from human genome, the
lengths of which are 1,000 bp each. Here we compare the mutation rates
across each type of region.

## Mutation rates in different genomic regions

The number of single mutations in each type of region:

| position |   count |   length | mutation\_rate |
| :------- | ------: | -------: | :------------- |
| enhancer |   20636 | 11759211 | 1.75e-05       |
| exon     |  126086 | 80790540 | 1.56e-05       |
| others   | 7363597 |       NA | NA             |
| promoter |   82772 | 46589695 | 1.78e-05       |
| random   |  237258 | 98325390 | 2.41e-05       |

## Mutation rates for different promoter length

Promoter regions are usually defined as the parts of genome adjacent to
each TSS. The exact length to the left and right of the TSS varies from
study to study. We want to investigate if the definition of promoter
will affect the mutation rate estimation, so we’ve tried different
length and computed the corresponding rate.

| promoter      |  count |    length | mutation\_rate |
| :------------ | -----: | --------: | :------------- |
| l1000\_r1000  |  82772 |  46589695 | 1.78e-05       |
| l1000\_r0     |  53791 |  29235415 | 1.84e-05       |
| l500\_r500    |  43214 |  21491714 | 2.01e-05       |
| l5000\_r1000  | 267209 | 153418487 | 1.74e-05       |
| l10000\_r1000 | 498132 | 275723250 | 1.81e-05       |

## Exon mutations

<!-- The transcripts whose exon harbor the greatest number of mutations: -->

The genes whose exons harbor the greatest number of
mutations:

| gene    | count | exon\_length | expected\_count |  var\_count | log.p\_value |
| :------ | ----: | -----------: | --------------: | ----------: | -----------: |
| TTN     |   402 |       105931 |     203.0683020 | 203.0594005 |   34.2297563 |
| MUC16   |   291 |        43900 |      38.8441462 |  38.8425092 |  146.5584208 |
| MUC19   |   145 |        25005 |     131.7101104 | 131.7046326 |    0.8755063 |
| DNAH5   |   130 |        15652 |      30.8350697 |  30.8336971 |   39.5101051 |
| GRIN2B  |   127 |        30361 |      11.6242257 |  11.6237793 |   83.1847409 |
| PTPRT   |   119 |        12668 |      12.4234055 |  12.4229130 |   71.8795748 |
| PCLO    |   110 |        17115 |       8.0240822 |   8.0237973 |   82.1696779 |
| LRP1B   |   106 |        16622 |      40.8092444 |  40.8074838 |   16.8350106 |
| FUT9    |    91 |        12801 |       8.1078761 |   8.1075213 |   60.9016136 |
| CSMD1   |    87 |        14387 |      22.6533751 |  22.6523481 |   24.1368090 |
| KSR2    |    85 |        17028 |      19.3304485 |  19.3296754 |   27.4041883 |
| ANK3    |    83 |        18148 |      19.2918264 |  19.2910455 |   26.1755776 |
| FAT4    |    81 |        19409 |       6.1076164 |   6.1073309 |   59.7265232 |
| ADGRV1  |    79 |        19423 |      72.3115418 |  72.3085136 |    0.6372575 |
| NEBL    |    78 |        15907 |      12.1527640 |  12.1522714 |   35.6548679 |
| DSCAM   |    77 |        10920 |      25.7172283 |  25.7160963 |   15.5710596 |
| CACNA1E |    76 |        15006 |      18.1954366 |  18.1946028 |   23.3038982 |
| PKHD1L1 |    74 |        20038 |      26.6169491 |  26.6157322 |   13.4287314 |
| GRIN2A  |    73 |        14463 |       7.6237776 |   7.6234151 |   44.5157142 |
| ZNF208  |    71 |        10832 |       0.9534276 |   0.9533854 |  103.8085164 |

## Promoter mutations

<!-- The transcripts whose promoters harbor the greatest number of mutations: -->

The genes whose promoters harbor the greatest number of
mutations:

| gene     | count | promoter\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| HLA-DRB5 |    87 |             2001 |        3.198050 |   3.197935 |    89.771631 |
| HLA-DRB1 |    82 |             2001 |        3.217392 |   3.217276 |    82.441833 |
| PRAMEF1  |    56 |             4002 |        6.486443 |   6.486214 |    32.144180 |
| TOLLIP   |    56 |             2001 |        3.859760 |   3.859618 |    43.650350 |
| FANK1    |    55 |             2001 |        3.765051 |   3.764909 |    43.041234 |
| CLEC5A   |    38 |             4002 |        6.021637 |   6.021437 |    17.632050 |
| UGT3A1   |    37 |             4002 |        6.508988 |   6.508752 |    15.784424 |
| TERT     |    35 |             2001 |        3.994264 |   3.994113 |    20.647739 |
| BAGE     |    33 |             2001 |        3.554660 |   3.554530 |    20.258326 |
| BAGE2    |    33 |             2001 |        3.563438 |   3.563308 |    20.226667 |
| BAGE3    |    33 |             2001 |        3.563438 |   3.563308 |    20.226667 |
| BAGE4    |    33 |             2001 |        3.563438 |   3.563308 |    20.226667 |
| BAGE5    |    33 |             2001 |        3.563438 |   3.563308 |    20.226667 |
| MYH2     |    32 |             4002 |        5.742037 |   5.741838 |    13.541327 |
| CWH43    |    31 |             4002 |        7.122240 |   7.121975 |    10.468756 |
| NLRP11   |    31 |             4002 |        6.502088 |   6.501862 |    11.436387 |
| PLCZ1    |    31 |             4002 |        5.834469 |   5.834278 |    12.616074 |
| TPTE     |    31 |             2001 |        3.876607 |   3.876462 |    17.300770 |
| HLA-A    |    30 |             2001 |        3.825555 |   3.825408 |    16.547307 |
| KRT18    |    30 |             4002 |        7.876875 |   7.876562 |     8.828112 |

## Enhancer mutations

<!-- The transcripts whose enhancers harbor the greatest number of mutations: -->

The genes whose enhancers harbor the greatest number of
mutations:

| gene     | count | enhancer\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| CD86     |    62 |             7857 |        12.75901 |   12.75837 |    22.380567 |
| FPR1     |    58 |            11798 |        20.84897 |   20.84786 |    10.731581 |
| ST3GAL1  |    56 |            18914 |        33.16872 |   33.16695 |     3.729844 |
| HLA-DMA  |    55 |             8575 |        14.41119 |   14.41045 |    15.505565 |
| HLA-DMB  |    55 |             8564 |        14.40643 |   14.40568 |    15.511446 |
| HLA-DRB5 |    55 |             8373 |        14.16163 |   14.16090 |    15.817001 |
| IL1A     |    55 |            18405 |        31.90138 |   31.89970 |     3.894158 |
| BRD2     |    54 |             8752 |        14.62609 |   14.62533 |    14.665155 |
| HLA-DQB1 |    54 |             9358 |        15.74156 |   15.74074 |    13.413991 |
| HLA-DRB1 |    54 |             9202 |        15.51845 |   15.51765 |    13.654271 |
| IL1RN    |    54 |            18697 |        32.46350 |   32.46179 |     3.472948 |
| IL36G    |    54 |            18697 |        32.46350 |   32.46179 |     3.472948 |
| FPR2     |    53 |            10836 |        19.23599 |   19.23496 |     9.737972 |
| HLA-DOA  |    53 |             8466 |        14.21495 |   14.21421 |    14.577206 |
| HLA-DOB  |    53 |             8438 |        14.11435 |   14.11362 |    14.698061 |
| HLA-DPB1 |    53 |             9300 |        15.85636 |   15.85553 |    12.756812 |
| SIGLEC10 |    53 |            11135 |        20.26162 |   20.26052 |     8.975130 |
| TAS2R38  |    53 |             6830 |        10.77307 |   10.77255 |    19.499454 |
| IL36B    |    52 |            18202 |        31.59279 |   31.59113 |     3.269738 |
| ZNF350   |    52 |            10382 |        18.51218 |   18.51119 |     9.854186 |

## Regulatory mutations

<!-- The transcripts whose regulatory regions have the most significant excess of mutations: -->

The genes whose regulatory regions have the most significant excess of
mutations:

| gene     | count | promoter\_length | enhancer\_length | total\_length | expected\_count | log.p\_value |
| :------- | ----: | ---------------: | ---------------: | ------------: | --------------: | -----------: |
| HLA-DRB5 |   142 |             2001 |             8373 |         10374 |       17.359684 |     76.89874 |
| HLA-DRB1 |   136 |             2001 |             9202 |         11203 |       18.735846 |     67.55288 |
| FANK1    |    55 |             2001 |               NA |          2001 |        3.765051 |     43.04123 |
| CLEC5A   |    88 |             4002 |             6339 |         10341 |       16.054778 |     35.06146 |
| CD86     |    90 |             4002 |             7857 |         11859 |       18.799167 |     31.56365 |
| PRAMEF1  |    56 |             4002 |              172 |          4174 |        6.878431 |     30.88402 |
| TOLLIP   |    57 |             2001 |             2651 |          4652 |        9.503698 |     24.91783 |
| TAS2R38  |    67 |             2001 |             6830 |          8831 |       13.846112 |     24.00780 |
| BAGE     |    38 |             2001 |              799 |          2800 |        4.519643 |     21.73411 |
| BAGE2    |    38 |             2001 |              799 |          2800 |        4.528421 |     21.70579 |
| BAGE3    |    33 |             2001 |               NA |          2001 |        3.563438 |     20.22667 |
| BAGE4    |    33 |             2001 |               NA |          2001 |        3.563438 |     20.22667 |
| BAGE5    |    33 |             2001 |               NA |          2001 |        3.563438 |     20.22667 |
| HLA-DRA  |    55 |             2001 |             5092 |          7093 |       11.450270 |     19.74275 |
| TERT     |    39 |             2001 |             1060 |          3061 |        5.849224 |     18.86468 |
| GIMAP4   |    50 |             2001 |             4149 |          6150 |       10.094181 |     18.56808 |
| DEFB121  |    38 |             2001 |             1517 |          3518 |        5.827790 |     18.09039 |
| TPTE     |    31 |             2001 |               NA |          2001 |        3.876607 |     17.30077 |
| HLA-DQB1 |    67 |             2001 |             9358 |         11359 |       18.915924 |     17.08893 |
| OR5AK2   |    29 |             2001 |              395 |          2396 |        3.429968 |     16.86024 |
