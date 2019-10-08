# ICGC Cancer Mutations
A deep dive into ICGC cancer mutations

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

Simple somatic mutation data submitted by member projects are annotated using the variant annotation tool [SnpEff](http://snpeff.sourceforge.net/SnpEff.html). SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes) using Sequence Ontology terms by default. It has a web-based user interface via [Galaxy](https://usegalaxy.org) project.

SnpEff can annotate SNPs, MNPs, insertions and deletions. Support for mixed variants and structural variants is available (although sometimes limited). It also provides a simple assessment of the putative impact of the variant (e.g. HIGH, MODERATE or LOW impact).

A typical SnpEff use case would be:

* Input: The inputs are predicted variants (SNPs, insertions, deletions and MNPs). The input file is usually obtained as a result of a sequencing experiment, and it is usually in variant call format (VCF).
* Output: SnpEff analyzes the input variants. It annotates the variants and calculates the effects they produce on known genes (e.g. amino acid changes).

#### VCF files
Variant Call Format (VCF) is the recomended format for input files. This is the format used by the "1000 Genomes Project", and is currently considered the de-facto standard for genomic variants. It is also the default format used in SnpEff.

In a nutshell, VCF format is tab-separated text file having the following columns:

1. Chromosome name
2. Position
3. Variant's ID
4. Reference genome
5. Alternative (i.e. variant)
6. Quality score
7. Filter (whether or not the variant passed quality filters)
8. INFO : Generic information about this variant. SnpEff adds annotation information in this column.

Here is an example of a few lines in a VCF file:

```
#CHROM POS     ID        REF    ALT     QUAL FILTER INFO                    
20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2
20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017   
```

Note that the first line is header information. Header lines start with '#'

#### Running SnpEff
SnpEff is very easy to install

```
# Download using wget
$ wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip

# If you prefer to use 'curl' instead of 'wget', you can type:
#     curl -L http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip > snpEff_latest_core.zip

# Install
$ unzip snpEff_latest_core.zip 
```

After installation, the basic command to run an annotation is 

```
$ java -Xmx4g -jar snpEff.jar GRCh37.75 path/to/input/input.vcf > path/to/output/output.vcf
```

You can run SnpEff in "the Cloud" exactly the same way as running it on your local computer. You should not have any problems at all.
Here is an example of installing it and running it on an Amazon EC2 instance (virtual machine):

```
$ ssh -i ./aws_amazon/pcingola_aws.pem ec2-user@ec2-54-234-14-244.compute-1.amazonaws.com

       __|  __|_  )
       _|  (     /   Amazon Linux AMI
      ___|\___|___|


[ec2-user@ip-10-2-202-163 ~]$ wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
[ec2-user@ip-10-2-202-163 ~]$ unzip snpEff_latest_core.zip
[ec2-user@ip-10-2-202-163 ~]$ cd snpEff/
[ec2-user@ip-10-2-202-163 snpEff]$ java -jar snpEff.jar download -v hg19
00:00:00.000    Downloading database for 'hg19'
...
00:00:36.340    Done
[ec2-user@ip-10-2-202-163 snpEff]$ java -Xmx4G -jar snpEff.jar dump -v hg19 > /dev/null
00:00:00.000    Reading database for genome 'hg19' (this might take a while)
00:00:20.688    done
00:00:20.688    Building interval forest
00:00:33.110    Done.
```
More options can be found at the [SnpEff manual](http://snpeff.sourceforge.net/SnpEff_manual.html)

### Mutation function impact prediction
To assess the functional impact of missense somatic mutations on protein coding genes, ICGC uses FatHMM to compute functional impact scores and assign impact categories. ICGC DCC defines the following four categories: High, Medium, Low and Unknown. This is not very related to our study.

### Simple somatic mutation frequencies
In order to assess how often a particular simple somatic mutation (SSM) occurs in cancer patients, we establishes specific criteria to determine whether SSMs reported from different donors are considered to be the same. When all of these fields from two reported SSMs: chromosome, chromosome start, chromosome end, mutation type, DNA change, genome assembly version are identical, they are merged into one mutation entity. With this, we build a non-redundant collection of all simple somatic mutations. Each mutation entry in the collection is assigned a stable identifier for persisted referencing portal wide and across releases. Information for linking these non-redundant SSMs to the reported ones is also kept, enabling mutation counts across donors.

## Explorative data analysis
Each cancer project has some basic numbers reported on the [Cancer Projects](https://dcc.icgc.org/projects) page, namely the number of donors, the number of mutations, and the number of somatic mutations in donor's exomes. We reproduced these numbers. The details can be found in the [icgc-eda](https://github.com/jz132/cancer-mutations/tree/master/icgc-eda) directory.

## Construct enhancer regions and connect them to target genes
Since we are more interested in the noncoding somatic mutations in cancer, an important step in our study is to focus on mutations in noncoding regions that are likely to have functional impact on gene expression. In particular, we want to investigate the mutations in promoter and enhancer regions.

While promoters can be easily defined as the genomic regions close to TSS, the definition of enhancers are quite ambiguous and versatile. Moreover, since enhancers can be far away from the gene they regulate, it is crucial to identify the connection between them. Here we adapt two recently published work to construct enhancers and link them to TSS, [An atlas of active enhancers across human cell types and tissues](https://www.nature.com/articles/nature12787) and [Mapping and analysis of chromatin state dynamics in nine human cell types](https://www.nature.com/articles/nature09906).

The details can be found in the [pelinks](https://github.com/jz132/cancer-mutations/tree/master/pelinks) directory.

