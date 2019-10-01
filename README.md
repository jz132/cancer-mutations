# ICGC Cancer Mutations
a deep dive into ICGC cancer mutations

## About the ICGC data portal
Most of the information about the [ICGC data portal](https://dcc.icgc.org) can be found at https://docs.icgc.org/portal/about/. Here we only list the information related to simple somatic mutations data. 

### Gene model
The canonical gene data (release 75) is available from Ensembl FTP site at ftp://ftp.ensembl.org/pub/release-75/mysql/homo_sapiens_core_75_37. The files related to ICGC are:

* gene.txt.gz
* transcript.txt.gz
* translation.txt.gz
* exon.txt.gz
* xref.txt.gz

However, since the headers are not included in these files, it is not clear to me how these data are integrated together. These efforts generate the information in the Genes tab on the website. 

Anyway, I don't think the gene annotation is as important as the mutation annotation for our study.

### Mutation consequence annotation
ICGC uses [Sequence Ontology's](www.sequenceontology.org) controlled vocabulary regarding mutation-induced changes. 

Simple somatic mutation data submitted by member projects are annotated using the variant annotation tool [SnpEff](http://snpeff.sourceforge.net/SnpEff.html). SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes).

A typical SnpEff use case would be:

* Input: The inputs are predicted variants (SNPs, insertions, deletions and MNPs). The input file is usually obtained as a result of a sequencing experiment, and it is usually in variant call format (VCF).
* Output: SnpEff analyzes the input variants. It annotates the variants and calculates the effects they produce on known genes (e.g. amino acid changes).

