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

| gene    | count | exon\_length | expected\_count | var\_count | log.p\_value |
| :------ | ----: | -----------: | --------------: | ---------: | -----------: |
| ZNF460  |     4 |         6882 |       0.0092932 |  0.0092931 |     9.510787 |
| DIO3    |     2 |         2103 |       0.0011098 |  0.0011098 |     6.210889 |
| TRIM71  |     2 |         8689 |       0.0054172 |  0.0054172 |     4.835043 |
| CPEB2   |     3 |         7182 |       0.0473501 |  0.0473500 |     4.767593 |
| PRSS23  |     2 |         7612 |       0.0062374 |  0.0062373 |     4.712835 |
| B3GAT1  |     2 |         3531 |       0.0064964 |  0.0064964 |     4.677559 |
| ZNF703  |     2 |         3351 |       0.0070884 |  0.0070884 |     4.601983 |
| NECTIN3 |     2 |         9161 |       0.0077271 |  0.0077271 |     4.527231 |
| MEGF9   |     2 |         6216 |       0.0077677 |  0.0077676 |     4.522698 |
| ZNF281  |     2 |         9666 |       0.0086292 |  0.0086291 |     4.431589 |
| AWAT1   |     2 |         1406 |       0.0090118 |  0.0090118 |     4.394014 |
| NUDT10  |     2 |         2656 |       0.0097834 |  0.0097834 |     4.322879 |
| DPM3    |     2 |          905 |       0.0103376 |  0.0103375 |     4.275185 |
| JRK     |     2 |        10808 |       0.0104491 |  0.0104491 |     4.265897 |
| BCL2L1  |     2 |         4949 |       0.0106247 |  0.0106247 |     4.251471 |
| TRPC3   |     3 |         4404 |       0.0768017 |  0.0768014 |     4.147007 |
| TET2    |     2 |        18880 |       0.0125279 |  0.0125279 |     4.108896 |
| LILRA3  |     2 |         1379 |       0.0126587 |  0.0126587 |     4.099914 |
| KDM6A   |     3 |         5970 |       0.0818289 |  0.0818287 |     4.066030 |
| USPL1   |     2 |         7590 |       0.0148309 |  0.0148309 |     3.962986 |

## Promoter mutations

<!-- The transcripts whose promoters harbor the greatest number of mutations: -->

The genes whose promoters harbor the greatest number of
mutations:

| gene       | count | promoter\_length | expected\_count | var\_count | log.p\_value |
| :--------- | ----: | ---------------: | --------------: | ---------: | -----------: |
| EFCAB3     |     3 |             4002 |       0.0375239 |  0.0375238 |     5.067439 |
| TRAPPC3L   |     2 |             2001 |       0.0183548 |  0.0183547 |     3.778842 |
| ZP4        |     2 |             2001 |       0.0185261 |  0.0185260 |     3.770823 |
| GP5        |     2 |             2001 |       0.0189266 |  0.0189265 |     3.752361 |
| LILRA3     |     2 |             2001 |       0.0191546 |  0.0191546 |     3.742025 |
| TAF13      |     2 |             2001 |       0.0192171 |  0.0192170 |     3.739214 |
| ALDOC      |     2 |             2001 |       0.0193976 |  0.0193976 |     3.731146 |
| CTSE       |     2 |             2001 |       0.0194712 |  0.0194711 |     3.727878 |
| PDGFD      |     2 |             2001 |       0.0199483 |  0.0199483 |     3.706988 |
| RPL22      |     2 |             2001 |       0.0202849 |  0.0202848 |     3.692554 |
| CELF2      |     2 |             2001 |       0.0202923 |  0.0202923 |     3.692237 |
| ST6GALNAC5 |     2 |             2001 |       0.0206028 |  0.0206028 |     3.679136 |
| RRP9       |     2 |             2001 |       0.0209410 |  0.0209410 |     3.665092 |
| PI4K2A     |     2 |             2001 |       0.0209484 |  0.0209484 |     3.664788 |
| PARP3      |     2 |             2001 |       0.0210738 |  0.0210737 |     3.659642 |
| SIRT7      |     2 |             2001 |       0.0212226 |  0.0212225 |     3.653573 |
| EPHA2      |     2 |             2001 |       0.0212546 |  0.0212545 |     3.652273 |
| SLC29A1    |     2 |             2001 |       0.0212820 |  0.0212819 |     3.651162 |
| GCSH       |     2 |             2001 |       0.0213751 |  0.0213751 |     3.647396 |
| HELT       |     2 |             2001 |       0.0214827 |  0.0214826 |     3.643067 |

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
