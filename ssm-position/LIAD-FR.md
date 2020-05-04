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
| l10000\_r1000 |  2796 | 275723250 | 2.03e-06       |

## Exon mutations

<!-- The transcripts whose exon harbor the greatest number of mutations: -->

The genes whose exons harbor the greatest number of
mutations:

| gene     | count | exon\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | -----------: | --------------: | ---------: | -----------: |
| ZNF460   |     4 |         6882 |       0.0717793 |  0.0717791 |     5.981123 |
| TRPC3    |     3 |         4404 |       0.0493490 |  0.0493489 |     4.714370 |
| KDM6A    |     3 |         5970 |       0.0618802 |  0.0618801 |     4.423620 |
| DPM3     |     2 |          905 |       0.0103783 |  0.0103783 |     4.271780 |
| CPEB2    |     3 |         7182 |       0.0783169 |  0.0783167 |     4.122044 |
| PGLYRP3  |     2 |         1180 |       0.0126851 |  0.0126851 |     4.098112 |
| AWAT1    |     2 |         1406 |       0.0143328 |  0.0143327 |     3.992516 |
| LILRA3   |     2 |         1379 |       0.0150442 |  0.0150441 |     3.950646 |
| DIO3     |     2 |         2103 |       0.0251124 |  0.0251123 |     3.508518 |
| IL18RAP  |     2 |         2683 |       0.0272223 |  0.0272222 |     3.439053 |
| KRT72    |     2 |         2483 |       0.0275262 |  0.0275261 |     3.429498 |
| NUDT10   |     2 |         2656 |       0.0308031 |  0.0308030 |     3.332747 |
| GMNC     |     2 |         3629 |       0.0364431 |  0.0364431 |     3.188334 |
| B3GAT1   |     2 |         3531 |       0.0398995 |  0.0398993 |     3.110629 |
| ZNF703   |     2 |         3351 |       0.0427163 |  0.0427162 |     3.052188 |
| C15orf41 |     2 |         4554 |       0.0460823 |  0.0460822 |     2.987279 |
| GAB3     |     2 |         4760 |       0.0486181 |  0.0486180 |     2.941481 |
| DMRT2    |     2 |         4645 |       0.0507013 |  0.0507011 |     2.905640 |
| BCL2L1   |     2 |         4949 |       0.0523391 |  0.0523389 |     2.878499 |
| SLC4A3   |     2 |         4525 |       0.0527090 |  0.0527089 |     2.872487 |

## Promoter mutations

<!-- The transcripts whose promoters harbor the greatest number of mutations: -->

The genes whose promoters harbor the greatest number of
mutations:

| gene     | count | promoter\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| EFCAB3   |     3 |             4002 |       0.0288322 |  0.0288322 |     5.407901 |
| TRAPPC3L |     2 |             2001 |       0.0142523 |  0.0142522 |     3.997387 |
| GP5      |     2 |             2001 |       0.0143733 |  0.0143733 |     3.990075 |
| TAF13    |     2 |             2001 |       0.0144971 |  0.0144971 |     3.982664 |
| ALDOC    |     2 |             2001 |       0.0146061 |  0.0146060 |     3.976190 |
| CTSE     |     2 |             2001 |       0.0147006 |  0.0147006 |     3.970613 |
| PDGFD    |     2 |             2001 |       0.0148392 |  0.0148392 |     3.962504 |
| RPL22    |     2 |             2001 |       0.0149862 |  0.0149862 |     3.953983 |
| CELF2    |     2 |             2001 |       0.0149987 |  0.0149987 |     3.953260 |
| RRP9     |     2 |             2001 |       0.0153536 |  0.0153536 |     3.933050 |
| PI4K2A   |     2 |             2001 |       0.0153980 |  0.0153980 |     3.930556 |
| PARP3    |     2 |             2001 |       0.0154400 |  0.0154399 |     3.928204 |
| SIRT7    |     2 |             2001 |       0.0154973 |  0.0154973 |     3.925002 |
| HELT     |     2 |             2001 |       0.0155742 |  0.0155742 |     3.920723 |
| TDRD9    |     2 |             2001 |       0.0156300 |  0.0156300 |     3.917634 |
| C17orf77 |     2 |             4002 |       0.0285958 |  0.0285958 |     3.396695 |
| MROH9    |     2 |             4002 |       0.0288234 |  0.0288233 |     3.389876 |
| CFP      |     2 |             4002 |       0.0295368 |  0.0295368 |     3.368844 |
| KCNK2    |     2 |             4002 |       0.0302886 |  0.0302885 |     3.347231 |
| TECR     |     2 |             4002 |       0.0304257 |  0.0304256 |     3.343347 |

## Enhancer mutations

<!-- The transcripts whose enhancers harbor the greatest number of mutations: -->

The genes whose enhancers harbor the greatest number of
mutations:

| gene     | count | enhancer\_length | expected\_count | var\_count | log.p\_value |
| :------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| CACNA2D1 |     2 |              922 |       0.0099831 |  0.0099831 |     4.305386 |
| NR4A2    |     2 |             2358 |       0.0286783 |  0.0286782 |     3.394217 |
| STX16    |     2 |             3376 |       0.0384574 |  0.0384572 |     3.142188 |
| ZNF831   |     2 |             5960 |       0.0660301 |  0.0660299 |     2.680611 |
| SNAPC1   |     1 |              234 |       0.0025863 |  0.0025863 |     2.587887 |
| RNF149   |     2 |             7118 |       0.0792278 |  0.0792275 |     2.526137 |
| TRPC5    |     1 |              287 |       0.0030512 |  0.0030511 |     2.516198 |
| ZCCHC16  |     1 |              287 |       0.0030512 |  0.0030511 |     2.516198 |
| KYNU     |     2 |             8048 |       0.0826284 |  0.0826281 |     2.490612 |
| MOG      |     1 |              344 |       0.0040303 |  0.0040303 |     2.395535 |
| OR2H1    |     1 |              480 |       0.0053784 |  0.0053783 |     2.270518 |
| RPAP3    |     1 |              541 |       0.0056543 |  0.0056543 |     2.248850 |
| MAP4K4   |     2 |            10046 |       0.1118534 |  0.1118530 |     2.235965 |
| PDE4DIP  |     2 |            10498 |       0.1125707 |  0.1125703 |     2.230618 |
| METTL22  |     1 |              611 |       0.0072344 |  0.0072344 |     2.142168 |
| TBC1D8   |     1 |              744 |       0.0079795 |  0.0079794 |     2.099758 |
| ZMYND8   |     1 |              699 |       0.0082470 |  0.0082469 |     2.085495 |
| RASAL1   |     1 |              672 |       0.0083282 |  0.0083282 |     2.081256 |
| AVP      |     1 |              915 |       0.0083748 |  0.0083748 |     2.078844 |
| OXT      |     1 |              915 |       0.0083748 |  0.0083748 |     2.078844 |

## Regulatory mutations

<!-- The transcripts whose regulatory regions have the most significant excess of mutations: -->

The genes whose regulatory regions have the most significant excess of
mutations:

| gene     | count | promoter\_length | enhancer\_length | total\_length | expected\_count | log.p\_value |
| :------- | ----: | ---------------: | ---------------: | ------------: | --------------: | -----------: |
| EFCAB3   |     3 |             4002 |             1121 |          5123 |       0.0421711 |     4.916827 |
| TRAPPC3L |     2 |             2001 |               NA |          2001 |       0.0142523 |     3.997387 |
| TAF13    |     2 |             2001 |               NA |          2001 |       0.0144971 |     3.982664 |
| ALDOC    |     2 |             2001 |               NA |          2001 |       0.0146061 |     3.976190 |
| CTSE     |     2 |             2001 |               NA |          2001 |       0.0147006 |     3.970613 |
| RPL22    |     2 |             2001 |               NA |          2001 |       0.0149862 |     3.953983 |
| RRP9     |     2 |             2001 |               NA |          2001 |       0.0153536 |     3.933050 |
| PARP3    |     2 |             2001 |               NA |          2001 |       0.0154400 |     3.928204 |
| HELT     |     2 |             2001 |               NA |          2001 |       0.0155742 |     3.920723 |
| TDRD9    |     2 |             2001 |               NA |          2001 |       0.0156300 |     3.917634 |
| GP5      |     2 |             2001 |              135 |          2136 |       0.0160204 |     3.896318 |
| MOG      |     2 |             2001 |              344 |          2345 |       0.0182599 |     3.783318 |
| RPAP3    |     2 |             2001 |              541 |          2542 |       0.0204832 |     3.684160 |
| PDGFD    |     2 |             2001 |              747 |          2748 |       0.0221015 |     3.618580 |
| CACNA2D1 |     2 |             2001 |              922 |          2923 |       0.0261949 |     3.472173 |
| C17orf77 |     2 |             4002 |               NA |          4002 |       0.0285958 |     3.396695 |
| MROH9    |     2 |             4002 |               NA |          4002 |       0.0288234 |     3.389876 |
| TECR     |     2 |             4002 |               NA |          4002 |       0.0304257 |     3.343347 |
| SNAP47   |     2 |             4002 |               NA |          4002 |       0.0310497 |     3.325893 |
| AMER3    |     2 |             4002 |               NA |          4002 |       0.0313887 |     3.316560 |
