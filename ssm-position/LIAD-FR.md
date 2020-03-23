LIAD-FR Simple Somatic Mutations
================

There are 5 donors with WGS data in the study. The total number of
single mutations is 40649. Here, we use enhancers from the FANTOM
project, exons from RefSeq, and promoters inferred from RefSeq TSS. The
three types of genomic region are defined so that there is no overlap
among them. We also take 100,000 random intervals from human genome, the
lengths of which are 1,000 bp each. Here we compare the mutation rates
across each type of region.

## Mutation rates in different genomic regions

The number of single mutations in each type of region:

| position | count |   length | mutation\_rate |
| :------- | ----: | -------: | :------------- |
| enhancer |   132 | 11759211 | 2.25e-06       |
| exon     |   853 | 80790540 | 2.11e-06       |
| others   | 39179 |       NA | NA             |
| promoter |   485 | 46589695 | 2.08e-06       |
| random   |  1289 | 98325390 | 2.62e-06       |

## Mutation rates for different promoter length

Promoter regions are usually defined as the parts of genome adjacent to
each TSS. The exact length to the left and right of the TSS varies from
study to study. We want to investigate if the definition of promoter
will affect the mutation rate estimation, so we’ve tried different
length and computed the corresponding rate.

| promoter      | count |    length | mutation\_rate |
| :------------ | ----: | --------: | :------------- |
| l1000\_r1000  |   485 |  46589695 | 2.08e-06       |
| l1000\_r0     |   296 |  29235415 | 2.02e-06       |
| l500\_r500    |   213 |  21491714 | 1.98e-06       |
| l5000\_r1000  |  1542 | 153418487 | 2.01e-06       |
| l10000\_r1000 |  2797 | 275723250 | 2.03e-06       |

## Exon mutations

<!-- The transcripts whose exon harbor the greatest number of mutations: -->

The genes whose exons harbor the greatest number of
mutations:

| gene     | count | exon\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | -----------: | --------------: | ---------: | -----------: |
| ZNF460   |     4 |         6882 |       0.0092932 |  0.0092931 |     9.510787 |
| CPEB2    |     3 |         7182 |       0.0473501 |  0.0473500 |     4.767593 |
| KDM6A    |     3 |         5970 |       0.0818289 |  0.0818287 |     4.066030 |
| TRPC3    |     3 |         4404 |       0.0768017 |  0.0768014 |     4.147007 |
| ACACA    |     2 |        10164 |       0.1684016 |  0.1684011 |     1.896750 |
| ATP8A2   |     2 |         9740 |       0.0802990 |  0.0802988 |     2.514780 |
| AWAT1    |     2 |         1406 |       0.0090118 |  0.0090118 |     4.394014 |
| B3GAT1   |     2 |         3531 |       0.0064964 |  0.0064964 |     4.677559 |
| BCL2L1   |     2 |         4949 |       0.0106247 |  0.0106247 |     4.251471 |
| BIRC6    |     2 |        15776 |       0.2661807 |  0.2661800 |     1.526879 |
| C15orf41 |     2 |         4554 |       0.0253536 |  0.0253536 |     3.500282 |
| CDH7     |     2 |         5728 |       0.0194990 |  0.0194990 |     3.726645 |
| COL5A2   |     2 |         6984 |       0.1156270 |  0.1156267 |     2.208227 |
| DDHD1    |     2 |        12911 |       0.0324550 |  0.0324549 |     3.287851 |
| DIO3     |     2 |         2103 |       0.0011098 |  0.0011098 |     6.210889 |
| DMD      |     2 |        14424 |       0.2264825 |  0.2264818 |     1.655909 |
| DMRT2    |     2 |         4645 |       0.0243103 |  0.0243102 |     3.536481 |
| DPM3     |     2 |          905 |       0.0103376 |  0.0103375 |     4.275185 |
| DSP      |     2 |        10319 |       0.0876518 |  0.0876515 |     2.440793 |
| FREM2    |     2 |        16186 |       0.0442429 |  0.0442428 |     3.022128 |

## Promoter mutations

<!-- The transcripts whose promoters harbor the greatest number of mutations: -->

The genes whose promoters harbor the greatest number of
mutations:

| gene     | count | promoter\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| EFCAB3   |     3 |             4002 |       0.0375239 |  0.0375238 |     5.067439 |
| ACTR3    |     2 |             4002 |       0.0421951 |  0.0421950 |     3.062701 |
| ALDOC    |     2 |             2001 |       0.0193976 |  0.0193976 |     3.731146 |
| AMER3    |     2 |             4002 |       0.0429743 |  0.0429742 |     3.047032 |
| C17orf77 |     2 |             4002 |       0.0375135 |  0.0375134 |     3.163499 |
| CDH7     |     2 |             4002 |       0.0419444 |  0.0419443 |     3.067805 |
| CELF2    |     2 |             2001 |       0.0202923 |  0.0202923 |     3.692237 |
| CFP      |     2 |             4002 |       0.0390422 |  0.0390421 |     3.129247 |
| CTSE     |     2 |             2001 |       0.0194712 |  0.0194711 |     3.727878 |
| DIO3     |     2 |             2001 |       0.0243678 |  0.0243678 |     3.534444 |
| DPM3     |     2 |             4002 |       0.0397448 |  0.0397447 |     3.113958 |
| EPHA2    |     2 |             2001 |       0.0212546 |  0.0212545 |     3.652273 |
| GCSH     |     2 |             2001 |       0.0213751 |  0.0213751 |     3.647396 |
| GP5      |     2 |             2001 |       0.0189266 |  0.0189265 |     3.752361 |
| HELT     |     2 |             2001 |       0.0214827 |  0.0214826 |     3.643067 |
| KCNK2    |     2 |             4002 |       0.0405311 |  0.0405310 |     3.097169 |
| LARP6    |     2 |             4002 |       0.0436841 |  0.0436840 |     3.033009 |
| LCN8     |     2 |             4002 |       0.0408851 |  0.0408850 |     3.089718 |
| LILRA3   |     2 |             2001 |       0.0191546 |  0.0191546 |     3.742025 |
| MROH9    |     2 |             4002 |       0.0370999 |  0.0370998 |     3.173011 |

## Enhancer mutations

<!-- The transcripts whose enhancers harbor the greatest number of mutations: -->

The genes whose enhancers harbor the greatest number of
mutations:

| gene     | count | enhancer\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| CACNA2D1 |     2 |              922 |       0.0099831 |  0.0099831 |    4.3053857 |
| FNDC3B   |     2 |            15941 |       0.1700519 |  0.1700514 |    1.8887504 |
| KYNU     |     2 |             8048 |       0.0826284 |  0.0826281 |    2.4906121 |
| MAP4K4   |     2 |            10046 |       0.1118534 |  0.1118530 |    2.2359646 |
| NR4A2    |     2 |             2358 |       0.0286783 |  0.0286782 |    3.3942173 |
| PDE4DIP  |     2 |            10498 |       0.1125707 |  0.1125703 |    2.2306184 |
| PLD1     |     2 |            13819 |       0.1472300 |  0.1472295 |    2.0074018 |
| RNF149   |     2 |             7118 |       0.0792278 |  0.0792275 |    2.5261375 |
| STX16    |     2 |             3376 |       0.0384574 |  0.0384572 |    3.1421879 |
| ZNF831   |     2 |             5960 |       0.0660301 |  0.0660299 |    2.6806112 |
| ABCD1    |     1 |             3065 |       0.0336767 |  0.0336766 |    1.4799628 |
| ADAM10   |     1 |             4883 |       0.0523322 |  0.0523320 |    1.2925456 |
| ADRB2    |     1 |             2488 |       0.0271941 |  0.0271940 |    1.5714175 |
| ANKRD28  |     1 |             5076 |       0.0554392 |  0.0554390 |    1.2681657 |
| APLF     |     1 |             5777 |       0.0622765 |  0.0622763 |    1.2191288 |
| ARHGAP25 |     1 |             4903 |       0.0538326 |  0.0538324 |    1.2805917 |
| ARHGAP4  |     1 |             2994 |       0.0319138 |  0.0319137 |    1.5029327 |
| ARHGEF1  |     1 |             8903 |       0.1056572 |  0.1056569 |    0.9988419 |
| ARID4B   |     1 |            10779 |       0.1213862 |  0.1213858 |    0.9419227 |
| ARID5B   |     1 |             8373 |       0.0855320 |  0.0855318 |    1.0863118 |

## Regulatory mutations

<!-- The transcripts whose regulatory regions have the most significant excess of mutations: -->

The genes whose regulatory regions have the most significant excess of
mutations:

| gene       | count | promoter\_length | enhancer\_length | total\_length | expected\_count | log.p\_value |
| :--------- | ----: | ---------------: | ---------------: | ------------: | --------------: | -----------: |
| EFCAB3     |     3 |             4002 |             1121 |          5123 |       0.0508628 |     4.675497 |
| TRAPPC3L   |     2 |             2001 |               NA |          2001 |       0.0183548 |     3.778842 |
| TAF13      |     2 |             2001 |               NA |          2001 |       0.0192171 |     3.739214 |
| ALDOC      |     2 |             2001 |               NA |          2001 |       0.0193976 |     3.731146 |
| CTSE       |     2 |             2001 |               NA |          2001 |       0.0194712 |     3.727878 |
| RPL22      |     2 |             2001 |               NA |          2001 |       0.0202849 |     3.692554 |
| GP5        |     2 |             2001 |              135 |          2136 |       0.0205737 |     3.680357 |
| RRP9       |     2 |             2001 |               NA |          2001 |       0.0209410 |     3.665092 |
| PARP3      |     2 |             2001 |               NA |          2001 |       0.0210738 |     3.659642 |
| SLC29A1    |     2 |             2001 |               NA |          2001 |       0.0212820 |     3.651162 |
| GCSH       |     2 |             2001 |               NA |          2001 |       0.0213751 |     3.647396 |
| HELT       |     2 |             2001 |               NA |          2001 |       0.0214827 |     3.643067 |
| TDRD9      |     2 |             2001 |               NA |          2001 |       0.0215990 |     3.638410 |
| ZP4        |     2 |             2001 |              376 |          2377 |       0.0221917 |     3.615068 |
| MOG        |     2 |             2001 |              344 |          2345 |       0.0226301 |     3.598202 |
| DIO3       |     2 |             2001 |               NA |          2001 |       0.0243678 |     3.534444 |
| RPAP3      |     2 |             2001 |              541 |          2542 |       0.0255489 |     3.493674 |
| PDGFD      |     2 |             2001 |              747 |          2748 |       0.0272106 |     3.439422 |
| ST6GALNAC5 |     2 |             2001 |              862 |          2863 |       0.0298527 |     3.359694 |
| EPHA2      |     2 |             2001 |              880 |          2881 |       0.0318141 |     3.304989 |
