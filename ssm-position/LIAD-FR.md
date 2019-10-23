LIAD-FR Simple Somatic Mutations
================

There are 5 donors with WGS data in the study. The total number of
single mutations is 40649.

The number of single mutations in each type of region:

| position | count |   length | mutation\_rate |
| :------- | ----: | -------: | :------------- |
| enhancer |   132 | 11800463 | 2.24e-06       |
| exon     |   855 | 81028378 | 2.11e-06       |
| others   | 39177 |       NA | NA             |
| promoter |   485 | 46574079 | 2.08e-06       |

For promoter regions, we played with its length:

| promoter      | count |    length | mutation\_rate |
| :------------ | ----: | --------: | :------------- |
| l1000\_r1000  |   485 |  46574079 | 2.08e-06       |
| l1000\_r0     |   276 |  24143431 | 2.29e-06       |
| l500\_r500    |   213 |  21482947 | 1.98e-06       |
| l5000\_r1000  |  1454 | 144455017 | 2.01e-06       |
| l10000\_r1000 |  2637 | 258385125 | 2.04e-06       |

The top-mutated enhancers:

| enhancer                 | count | length |      rate |
| :----------------------- | ----: | -----: | --------: |
| chr1:172788027-172788361 |     2 |    334 | 0.0059880 |
| chr2:143809081-143809289 |     2 |    208 | 0.0096154 |
| chr20:57557787-57558211  |     2 |    424 | 0.0047170 |
| chr7:81858656-81858971   |     2 |    315 | 0.0063492 |
| chr8:90479154-90479488   |     2 |    334 | 0.0059880 |
| chr1:25361269-25361784   |     1 |    515 | 0.0019417 |
| chr1:26907360-26907709   |     1 |    349 | 0.0028653 |
| chr1:62072423-62072573   |     1 |    150 | 0.0066667 |
| chr1:117040346-117040781 |     1 |    435 | 0.0022989 |
| chr1:145012331-145013054 |     1 |    723 | 0.0013831 |

The genes whose enhancers harbor the greatest number of mutations:

| gene          | mut\_count | enh\_count |
| :------------ | ---------: | ---------: |
| ATP1A1        |          3 |         51 |
| LPP           |          3 |          7 |
| PDE4DIP       |          3 |         49 |
| CACNA2D1      |          2 |          3 |
| FNDC3B        |          2 |         44 |
| NM\_001134772 |          2 |         10 |
| NM\_001199241 |          2 |         21 |
| NM\_002662    |          2 |         32 |
| NM\_004834    |          2 |         33 |
| NM\_145687    |          2 |         33 |
