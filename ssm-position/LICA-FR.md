LICA-FR Simple Somatic Mutations
================

There are 49 donors with WGS data in the study. The total number of
single mutations is 1352766. Here, we use enhancers from the FANTOM
project, exons from RefSeq, and promoters inferred from RefSeq TSS. The
three types of genomic region are defined so that there is no overlap
among them. We also take 100,000 random intervals from human genome, the
lengths of which are 1,000 bp each. Here we compare the mutation rates
across each type of region.

## Mutation rates in different genomic regions

The number of single mutations in each type of region:

| position |   count |   length | mutation\_rate |
| :------- | ------: | -------: | :------------- |
| enhancer |    5340 | 11759211 | 9.27e-06       |
| exon     |   31457 | 80790540 | 7.95e-06       |
| others   | 1296032 |       NA | NA             |
| promoter |   19937 | 46589695 | 8.73e-06       |
| random   |   42263 | 98325390 | 8.77e-06       |

## Mutation rates for different promoter length

Promoter regions are usually defined as the parts of genome adjacent to
each TSS. The exact length to the left and right of the TSS varies from
study to study. We want to investigate if the definition of promoter
will affect the mutation rate estimation, so we’ve tried different
length and computed the corresponding rate.

| promoter      |  count |    length | mutation\_rate |
| :------------ | -----: | --------: | :------------- |
| l1000\_r1000  |  19937 |  46589695 | 8.73e-06       |
| l1000\_r0     |  12633 |  29235415 | 8.82e-06       |
| l500\_r500    |   9828 |  21491714 | 9.33e-06       |
| l5000\_r1000  |  61567 | 153418487 | 8.19e-06       |
| l10000\_r1000 | 110270 | 275723250 | 8.16e-06       |

## Exon mutations

<!-- The transcripts whose exon harbor the greatest number of mutations: -->

The genes whose exons harbor the greatest number of
mutations:

| gene     | count | exon\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | -----------: | --------------: | ---------: | -----------: |
| CTNNB1   |    14 |         4062 |       1.3489836 |  1.3489672 |     9.665504 |
| NBEA     |    20 |        13106 |       4.2622333 |  4.2621822 |     7.547142 |
| TRIM58   |    13 |         5162 |       1.9933841 |  1.9933558 |     6.699370 |
| GPRIN2   |    16 |         6901 |       3.1363463 |  3.1362960 |     6.652535 |
| KRT72    |    10 |         2483 |       1.1403332 |  1.1403151 |     6.437622 |
| TMEM121B |    16 |         8989 |       4.3396759 |  4.3396034 |     4.880613 |
| ZP4      |     7 |         1967 |       0.7466447 |  0.7466345 |     4.872836 |
| ZFP2     |     7 |         2433 |       0.7576664 |  0.7576579 |     4.832429 |
| DUSP22   |    10 |         4121 |       1.7630617 |  1.7630345 |     4.788151 |
| COL22A1  |    13 |         6424 |       3.0134392 |  3.0133893 |     4.771945 |
| CACNA1I  |    17 |         9933 |       4.9985461 |  4.9984619 |     4.703385 |
| OR5M1    |     5 |          949 |       0.3173090 |  0.3173051 |     4.686166 |
| ACKR1    |     9 |         3253 |       1.4580318 |  1.4580087 |     4.651718 |
| NPAP1    |    13 |         8054 |       3.1106460 |  3.1106021 |     4.631220 |
| MUC19    |    23 |        25005 |       8.4109393 |  8.4108372 |     4.611085 |
| SLC45A4  |    14 |         7992 |       3.7257921 |  3.7257325 |     4.440156 |
| VSTM2A   |     9 |         5068 |       1.6570617 |  1.6570414 |     4.228140 |
| FLRT2    |    12 |         8232 |       2.9424616 |  2.9424223 |     4.224692 |
| MUC16    |    33 |        43900 |      15.7769118 | 15.7767076 |     3.993117 |
| ZNF835   |     8 |         2776 |       1.3916485 |  1.3916255 |     3.990074 |

## Promoter mutations

<!-- The transcripts whose promoters harbor the greatest number of mutations: -->

The genes whose promoters harbor the greatest number of
mutations:

| gene     | count | promoter\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| TERT     |    18 |             2001 |       0.9047636 |  0.9047517 |    16.960507 |
| NPY4R    |     9 |             2001 |       0.6817949 |  0.6817870 |     7.322518 |
| CDH18    |     8 |             2001 |       0.5453041 |  0.5452989 |     6.922265 |
| OTOA     |     7 |             2001 |       0.5917625 |  0.5917565 |     5.521320 |
| ABCB5    |     8 |             4002 |       0.8882208 |  0.8882142 |     5.358476 |
| EPB41L3  |     9 |             4002 |       1.1933969 |  1.1933844 |     5.332485 |
| ZNF595   |     7 |             2001 |       0.6348671 |  0.6348603 |     5.323807 |
| ZNF718   |     7 |             2001 |       0.6353268 |  0.6353199 |     5.321780 |
| PNLIPRP1 |     8 |             4002 |       1.1040417 |  1.1040308 |     4.685084 |
| POLA1    |     8 |             4002 |       1.1350947 |  1.1350834 |     4.600537 |
| DMRT2    |     8 |             4002 |       1.1927707 |  1.1927585 |     4.450290 |
| LYPD6B   |     8 |             4002 |       1.2416401 |  1.2416267 |     4.329367 |
| FIGLA    |     6 |             2001 |       0.6605136 |  0.6605061 |     4.182418 |
| ZNF213   |     8 |             4002 |       1.3400835 |  1.3400684 |     4.101683 |
| TNFRSF19 |     6 |             2001 |       0.7335395 |  0.7335308 |     3.935993 |
| BCL2L15  |     5 |             2001 |       0.4893453 |  0.4893411 |     3.807144 |
| FCAR     |     5 |             2001 |       0.5046982 |  0.5046937 |     3.745550 |
| ROR2     |     6 |             2001 |       0.8057781 |  0.8057680 |     3.717739 |
| LCN8     |     8 |             4002 |       1.5319180 |  1.5318992 |     3.709598 |
| TFAP2D   |     5 |             2001 |       0.5296991 |  0.5296941 |     3.649494 |

## Enhancer mutations

<!-- The transcripts whose enhancers harbor the greatest number of mutations: -->

The genes whose enhancers harbor the greatest number of
mutations:

| gene    | count | enhancer\_length | expected\_count | var\_count | log.p\_value |
| :------ | ----: | ---------------: | --------------: | ---------: | -----------: |
| CCL2    |     9 |             3465 |       1.4943426 |  1.4943152 |     4.569551 |
| GABRB3  |     4 |              744 |       0.2968514 |  0.2968466 |     3.592672 |
| FRG1    |     7 |             3345 |       1.3164261 |  1.3164036 |     3.362017 |
| PAG1    |     8 |             4817 |       2.0008060 |  2.0007712 |     2.958808 |
| NLRP3   |    10 |             7075 |       3.0412497 |  3.0411952 |     2.914148 |
| RETN    |     4 |              899 |       0.4625293 |  0.4625197 |     2.879082 |
| CCL1    |     5 |             1960 |       0.8350302 |  0.8350150 |     2.769730 |
| PTGER4  |    11 |            10344 |       3.7461965 |  3.7461396 |     2.762055 |
| ST8SIA4 |     6 |             3508 |       1.2647283 |  1.2647087 |     2.710457 |
| CD300LD |    13 |            10511 |       4.9934967 |  4.9933987 |     2.699705 |
| C2CD2   |     8 |             4531 |       2.2140440 |  2.2140002 |     2.687027 |
| GALNT9  |     6 |             2192 |       1.2803944 |  1.2803656 |     2.684063 |
| C6orf27 |     3 |              409 |       0.2488728 |  0.2488674 |     2.670771 |
| PRKAA1  |     4 |             1583 |       0.5573586 |  0.5573505 |     2.587462 |
| XAB2    |     5 |             1910 |       0.9445051 |  0.9444862 |     2.540936 |
| EIF3A   |     5 |             1991 |       0.9667405 |  0.9667213 |     2.498257 |
| CCL8    |    11 |             9358 |       4.0705528 |  4.0704779 |     2.490336 |
| ASH2L   |     2 |              224 |       0.0826581 |  0.0826569 |     2.490309 |
| POSTN   |     3 |              738 |       0.2909791 |  0.2909743 |     2.480643 |
| GNAI3   |     7 |             3890 |       1.9021740 |  1.9021360 |     2.459992 |

## Regulatory mutations

<!-- The transcripts whose regulatory regions have the most significant excess of mutations: -->

The genes whose regulatory regions have the most significant excess of
mutations:

| gene     | count | promoter\_length | enhancer\_length | total\_length | expected\_count | log.p\_value |
| :------- | ----: | ---------------: | ---------------: | ------------: | --------------: | -----------: |
| TERT     |    18 |             2001 |             1060 |          3061 |       1.4792078 |    13.353180 |
| NPY4R    |     9 |             2001 |               NA |          2001 |       0.6817949 |     7.322518 |
| CDH18    |     8 |             2001 |              151 |          2152 |       0.6070814 |     6.573089 |
| CCL2     |    12 |             2001 |             3465 |          5466 |       2.0502997 |     5.755453 |
| ABCB5    |     8 |             4002 |               NA |          4002 |       0.8882208 |     5.358476 |
| ZNF595   |     7 |             2001 |               NA |          2001 |       0.6348671 |     5.323807 |
| ZNF718   |     7 |             2001 |               NA |          2001 |       0.6353268 |     5.321780 |
| OTOA     |     7 |             2001 |              117 |          2118 |       0.6357745 |     5.319807 |
| PNLIPRP1 |     9 |             4002 |              509 |          4511 |       1.3082371 |     5.017702 |
| POLA1    |     8 |             4002 |               NA |          4002 |       1.1350947 |     4.600537 |
| TNFRSF19 |     7 |             2001 |              368 |          2369 |       0.8572559 |     4.494396 |
| LYPD6B   |     8 |             4002 |               NA |          4002 |       1.2416401 |     4.329367 |
| DMRT2    |     8 |             4002 |              381 |          4383 |       1.3186418 |     4.149580 |
| ZNF213   |     8 |             4002 |               NA |          4002 |       1.3400835 |     4.101683 |
| EPB41L3  |    12 |             4002 |             4632 |          8634 |       3.0392463 |     4.094068 |
| POSTN    |     6 |             2001 |              738 |          2739 |       0.7090115 |     4.015609 |
| CCL1     |     8 |             2001 |             1960 |          3961 |       1.3991466 |     3.974249 |
| ST6GAL2  |     9 |             4002 |              659 |          4661 |       1.8451693 |     3.880067 |
| TTBK1    |     7 |             2001 |              755 |          2756 |       1.1011391 |     3.824593 |
| BCL2L15  |     5 |             2001 |               NA |          2001 |       0.4893453 |     3.807144 |
