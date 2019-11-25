LIAD-FR Simple Somatic Mutations
================

There are 5 donors with WGS data in the study. The total number of
single mutations is 40649. Here, we use enhancers from the FANTOM
project, exons from RefSeq, and promoters inferred from RefSeq tss. The
three types of genomic regions are defined so that there is no overlap
among them. We also take 100,000 random intervals from human genome.
Each interval is 1,000 bp long.

The number of single mutations in each type of region:

| position | count |   length | mutation\_rate |
| :------- | ----: | -------: | :------------- |
| enhancer |   132 | 11759209 | 2.25e-06       |
| exon     |   853 | 80790540 | 2.11e-06       |
| others   | 39179 |       NA | NA             |
| promoter |   485 | 46589695 | 2.08e-06       |
| random   |  1289 | 98325390 | 2.62e-06       |

For promoter regions, we played with its length:

| promoter      | count |    length | mutation\_rate |
| :------------ | ----: | --------: | :------------- |
| l1000\_r1000  |   485 |  46589695 | 2.08e-06       |
| l1000\_r0     |   276 |  24139178 | 2.29e-06       |
| l500\_r500    |   213 |  21491714 | 1.98e-06       |
| l5000\_r1000  |  1454 | 144491937 | 2.01e-06       |
| l10000\_r1000 |  2638 | 258443502 | 2.04e-06       |

The top-mutated
enhancers:

| chromosome |     start |       end | enhancer                 | count | length | mut\_rate |
| :--------- | --------: | --------: | :----------------------- | ----: | -----: | --------: |
| chr1       | 172788028 | 172788361 | chr1:172788027-172788361 |     2 |    333 | 0.0060060 |
| chr2       | 143809082 | 143809289 | chr2:143809081-143809289 |     2 |    207 | 0.0096618 |
| chr20      |  57557788 |  57558211 | chr20:57557787-57558211  |     2 |    423 | 0.0047281 |
| chr7       |  81858657 |  81858971 | chr7:81858656-81858971   |     2 |    314 | 0.0063694 |
| chr8       |  90479155 |  90479488 | chr8:90479154-90479488   |     2 |    333 | 0.0060060 |
| chr1       |  25361270 |  25361784 | chr1:25361269-25361784   |     1 |    514 | 0.0019455 |
| chr1       |  26907361 |  26907709 | chr1:26907360-26907709   |     1 |    348 | 0.0028736 |
| chr1       |  62072424 |  62072573 | chr1:62072423-62072573   |     1 |    149 | 0.0067114 |
| chr1       | 117040347 | 117040781 | chr1:117040346-117040781 |     1 |    434 | 0.0023041 |
| chr1       | 145012332 | 145013054 | chr1:145012331-145013054 |     1 |    722 | 0.0013850 |

The tss whose enhancers harbor the greatest number of mutations:

| tss           | mut\_count | e\_count |
| :------------ | ---------: | -------: |
| NM\_000722    |          2 |        3 |
| NM\_001001433 |          2 |       10 |
| NM\_001032998 |          2 |       21 |
| NM\_001130081 |          2 |       32 |
| NM\_001134772 |          2 |       10 |
| NM\_001134773 |          2 |       10 |
| NM\_001199241 |          2 |       21 |
| NM\_001204868 |          2 |       10 |
| NM\_001242559 |          2 |       33 |
| NM\_001242560 |          2 |       33 |

The genes whose enhancers harbor the greatest number of mutations:

| gene     | mut\_count | e\_count |
| :------- | ---------: | -------: |
| CACNA2D1 |          2 |        3 |
| FNDC3B   |          2 |       39 |
| KYNU     |          2 |       21 |
| MAP4K4   |          2 |       33 |
| NR4A2    |          2 |        6 |
| PDE4DIP  |          2 |       26 |
| PLD1     |          2 |       32 |
| RNF149   |          2 |       21 |
| STX16    |          2 |       10 |
| ZNF831   |          2 |       17 |
