BRCA-EU Simple Somatic Mutations
================

There are 569 donors with WGS data in the study. The total number of
single mutations is 3594783. Here, we use enhancers from the FANTOM
project, exons from RefSeq, and promoters inferred from RefSeq tss. The
three types of genomic regions are defined so that there is no overlap
among them. We also take 100,000 random intervals from human genome.
Each interval is 1,000 bp long.

The number of single mutations in each type of region:

| position |   count |   length | mutation\_rate |
| :------- | ------: | -------: | :------------- |
| enhancer |   13648 | 11759209 | 2.04e-06       |
| exon     |   89381 | 80790540 | 1.94e-06       |
| others   | 3440711 |       NA | NA             |
| promoter |   51043 | 46589695 | 1.93e-06       |
| random   |  105777 | 98325390 | 1.89e-06       |

For promoter regions, we played with its length:

| promoter      |  count |    length | mutation\_rate |
| :------------ | -----: | --------: | :------------- |
| l1000\_r1000  |  51043 |  46589695 | 1.93e-06       |
| l1000\_r0     |  26497 |  24139178 | 1.93e-06       |
| l500\_r500    |  23992 |  21491714 | 1.96e-06       |
| l5000\_r1000  | 156717 | 144491937 | 1.91e-06       |
| l10000\_r1000 | 281208 | 258443502 | 1.91e-06       |

The top-mutated
enhancers:

| chromosome |     start |       end | enhancer                  | count | length | mut\_rate |
| :--------- | --------: | --------: | :------------------------ | ----: | -----: | --------: |
| chr1       | 177903072 | 177903284 | chr1:177903071-177903284  |    12 |    212 | 0.0566038 |
| chr15      |  65596703 |  65597188 | chr15:65596702-65597188   |     8 |    485 | 0.0164948 |
| chr10      |  30831189 |  30831972 | chr10:30831188-30831972   |     6 |    783 | 0.0076628 |
| chr11      | 129512960 | 129513316 | chr11:129512959-129513316 |     6 |    356 | 0.0168539 |
| chr12      |  69731954 |  69732266 | chr12:69731953-69732266   |     6 |    312 | 0.0192308 |
| chr12      |  70634013 |  70634396 | chr12:70634012-70634396   |     6 |    383 | 0.0156658 |
| chr17      |  38477237 |  38480096 | chr17:38477236-38480096   |     6 |   2859 | 0.0020986 |
| chr1       | 160677743 | 160679048 | chr1:160677742-160679048  |     5 |   1305 | 0.0038314 |
| chr10      |  30872539 |  30873163 | chr10:30872538-30873163   |     5 |    624 | 0.0080128 |
| chr10      |  38202648 |  38203041 | chr10:38202647-38203041   |     5 |    393 | 0.0127226 |

The tss whose enhancers harbor the greatest number of mutations:

| tss           | mut\_count | e\_count |
| :------------ | ---------: | -------: |
| NM\_005194    |         39 |       78 |
| NM\_005985    |         36 |       70 |
| NM\_001162505 |         35 |       67 |
| NM\_199129    |         35 |       67 |
| NM\_199203    |         35 |       67 |
| NM\_001731    |         33 |       52 |
| NM\_006290    |         32 |       54 |
| NM\_001778    |         31 |       33 |
| NM\_001160124 |         30 |       57 |
| NM\_001160125 |         30 |       57 |

The genes whose enhancers harbor the greatest number of mutations:

| gene    | mut\_count | e\_count |
| :------ | ---------: | -------: |
| CEBPB   |         39 |       78 |
| SNAI1   |         36 |       70 |
| TMEM189 |         35 |       67 |
| BTG1    |         33 |       52 |
| TNFAIP3 |         32 |       54 |
| CD48    |         31 |       33 |
| CD244   |         30 |       26 |
| KLF6    |         30 |       57 |
| SLC9A8  |         29 |       54 |
| FNDC3B  |         27 |       39 |
