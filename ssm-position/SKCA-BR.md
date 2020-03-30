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

| gene     | count | exon\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | -----------: | --------------: | ---------: | -----------: |
| MUC16    |   291 |        43900 |       71.784677 |  71.781637 |     83.19240 |
| TTN      |   402 |       105931 |      151.380223 | 151.374066 |     63.16541 |
| DNAH5    |   130 |        15652 |       23.037789 |  23.036784 |     52.61472 |
| HLA-DRB1 |    55 |         1468 |        2.567330 |   2.567211 |     51.67680 |
| PTPRT    |   119 |        12668 |       20.486491 |  20.485605 |     49.49759 |
| FUT9     |    91 |        12801 |       15.262553 |  15.261998 |     38.97075 |
| PAK5     |    65 |         4788 |        7.185567 |   7.185257 |     38.31705 |
| MUC19    |   145 |        25005 |       39.578032 |  39.576259 |     37.32643 |
| PCLO     |   110 |        17115 |       24.833783 |  24.832759 |     35.42171 |
| LRP1B    |   106 |        16622 |       23.261692 |  23.260705 |     35.19159 |
| FAM83B   |    59 |         6249 |        8.026848 |   8.026538 |     30.19768 |
| MGAM2    |    67 |         7920 |       11.268676 |  11.268197 |     28.90195 |
| MST1L    |    53 |         4683 |        6.786428 |   6.786158 |     28.44308 |
| ZNF208   |    71 |        10832 |       13.302663 |  13.302200 |     27.81888 |
| UGT3A1   |    55 |         5236 |        7.985758 |   7.985429 |     26.87785 |
| VSTM2A   |    51 |         5068 |        6.732631 |   6.732361 |     26.81713 |
| ST6GAL2  |    58 |         6926 |        9.561233 |   9.560857 |     25.57724 |
| TPTE     |    39 |         2789 |        3.864567 |   3.864400 |     25.04701 |
| ANKRD30A |    46 |         4441 |        5.813649 |   5.813414 |     25.04359 |
| SGCZ     |    55 |         7215 |        8.840980 |   8.840645 |     24.81139 |

## Promoter mutations

<!-- The transcripts whose promoters harbor the greatest number of mutations: -->

The genes whose promoters harbor the greatest number of
mutations:

| gene     | count | promoter\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| HLA-DRB5 |    87 |             2001 |        3.198050 |   3.197935 |     89.77163 |
| HLA-DRB1 |    82 |             2001 |        3.217392 |   3.217276 |     82.44183 |
| TOLLIP   |    56 |             2001 |        3.859760 |   3.859618 |     43.65035 |
| FANK1    |    55 |             2001 |        3.765051 |   3.764909 |     43.04123 |
| PRAMEF1  |    56 |             4002 |        6.486443 |   6.486214 |     32.14418 |
| TERT     |    35 |             2001 |        3.994264 |   3.994113 |     20.64774 |
| BAGE     |    33 |             2001 |        3.554660 |   3.554530 |     20.25833 |
| BAGE2    |    33 |             2001 |        3.563438 |   3.563308 |     20.22667 |
| BAGE3    |    33 |             2001 |        3.563438 |   3.563308 |     20.22667 |
| BAGE4    |    33 |             2001 |        3.563438 |   3.563308 |     20.22667 |
| BAGE5    |    33 |             2001 |        3.563438 |   3.563308 |     20.22667 |
| OR5AK2   |    28 |             2001 |        2.738911 |   2.738823 |     18.37852 |
| CLEC5A   |    38 |             4002 |        6.021637 |   6.021437 |     17.63205 |
| TPTE     |    31 |             2001 |        3.876607 |   3.876462 |     17.30077 |
| HLA-A    |    30 |             2001 |        3.825555 |   3.825408 |     16.54731 |
| UGT3A1   |    37 |             4002 |        6.508988 |   6.508752 |     15.78442 |
| PRAMEF2  |    26 |             2001 |        3.031837 |   3.031736 |     15.34649 |
| PRDM9    |    28 |             2001 |        3.695050 |   3.694913 |     15.13662 |
| OR4A16   |    24 |             2001 |        2.831572 |   2.831482 |     14.12184 |
| MYH2     |    32 |             4002 |        5.742037 |   5.741838 |     13.54133 |

## Enhancer mutations

<!-- The transcripts whose enhancers harbor the greatest number of mutations: -->

The genes whose enhancers harbor the greatest number of
mutations:

| gene     | count | enhancer\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| HOOK3    |    32 |              715 |       0.6758327 |  0.6758164 |     41.14985 |
| PADI2    |    33 |              951 |       1.9769904 |  1.9768766 |     28.00316 |
| CD86     |    62 |             7857 |      12.7590055 | 12.7583720 |     22.38057 |
| HLA-DRA  |    50 |             5092 |       8.3587083 |  8.3582865 |     21.92883 |
| TAS2R38  |    53 |             6830 |      10.7730724 | 10.7725458 |     19.49945 |
| CLEC5A   |    50 |             6339 |      10.0331411 | 10.0326487 |     18.67392 |
| ZNF597   |    21 |              564 |       1.4119759 |  1.4118763 |     17.14646 |
| NLRC3    |    21 |              621 |       1.4778359 |  1.4777340 |     16.75791 |
| EAF2     |    45 |             6038 |       9.5209238 |  9.5204635 |     16.07204 |
| DEFB116  |    25 |             1478 |       2.5597142 |  2.5595685 |     16.05271 |
| HLA-DRB5 |    55 |             8373 |      14.1616342 | 14.1608964 |     15.81700 |
| ZNF75A   |    21 |              726 |       1.6939286 |  1.6938184 |     15.60254 |
| HLA-DMB  |    55 |             8564 |      14.4064256 | 14.4056800 |     15.51145 |
| HLA-DMA  |    55 |             8575 |      14.4111919 | 14.4104466 |     15.50557 |
| PADI4    |    33 |             2735 |       5.5312666 |  5.5309518 |     14.75102 |
| HLA-DOB  |    53 |             8438 |      14.1143525 | 14.1136226 |     14.69806 |
| INHBA    |    48 |             7177 |      11.7691787 | 11.7685698 |     14.69093 |
| BRD2     |    54 |             8752 |      14.6260851 | 14.6253324 |     14.66515 |
| HLA-DOA  |    53 |             8466 |      14.2149454 | 14.2142101 |     14.57721 |
| DEFB121  |    23 |             1517 |       2.4915702 |  2.4914540 |     14.32831 |

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
