library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
library(fuzzyjoin)

## function reverseMerge ##
## combine the counts of a feature and its reverse complement.
## count_mat is a data frame of 4^k columns representing all k-mer features.
reverseMerge <- function(count_mat){
  features <- colnames(count_mat)
  list_drop <- NULL
  for(feature in features){
    reverse_feature <- as.character(reverseComplement(DNAString(feature)))
    if(feature < reverse_feature)
      count_mat[,feature] <- count_mat[,feature] + count_mat[,reverse_feature]
    else if(feature > reverse_feature)
      list_drop <- c(list_drop, which(features == feature))
    #if palindromic then feature == reverse_feature, and we count only once
  }
  count_mat <- count_mat[, -list_drop]
  return(count_mat)
}

# file paths and file names
icgc.data.path <- "/Users/jz132/Desktop/Gordanlab/Data/ICGC"
genomic.interval.path <- "/Users/jz132/r_projects/cancer-mutations/pelinks"
refseq.data.path <- "/Users/jz132/r_projects/cancer-mutations/pelinks/RefSeq"
output.path <- "/Users/jz132/r_projects/cancer-mutations/ssm-background/"
filename <- "simple_somatic_mutation.open.LIRI-JP.tsv"

genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg19")

keep_cols <- c("icgc_mutation_id", 
               "icgc_donor_id", 
               "chromosome",
               "chromosome_start",
               "chromosome_end",
               "assembly_version",
               "mutation_type",
               "reference_genome_allele",
               "mutated_from_allele",
               "mutated_to_allele",
               "consequence_type",
               "gene_affected",
               "transcript_affected",
               "sequencing_strategy")

# import ICGC data
setwd(icgc.data.path)
data_icgc_try <- read_tsv(filename, col_names = T, n_max = 5)
col_indicator <- paste(ifelse(colnames(data_icgc_try) %in% keep_cols, "?", "-"), collapse = "")
data_icgc_raw <- read_tsv(filename, col_names = T, col_types = col_indicator) %>%
  mutate(chromosome = paste0("chr", chromosome))

all_types <- data_icgc_raw %>% distinct(consequence_type) %>% pull()
nonexome <- c("intron_variant", "intragenic_variant", "upstream_gene_variant", 
              "downstream_gene_variant", "intergenic_region")
exome <- dplyr::setdiff(all_types, c(nonexome, NA))

# focus on single base substitution in WGS only
data_icgc_wgs <- data_icgc_raw %>% 
  filter(sequencing_strategy == "WGS", mutation_type == "single base substitution")
num_donors <- length(data_icgc_wgs %>% distinct(icgc_donor_id) %>% pull())
num_mutations <- length(data_icgc_wgs %>% distinct(icgc_mutation_id, icgc_donor_id) %>% pull())
if(!all(data_icgc_wgs$mutated_from_allele == data_icgc_wgs$reference_genome_allele)){
  warning("reference genome allele is not the same as mutated from allele")
}


# import genomic coordinates of exons
setwd(genomic.interval.path)
data_exons_fantom <- read_delim("all_exons_fantom.bed", delim = "\t", 
                                col_names = c("chromosome", "start", "end")) %>%
  arrange(chromosome, start, end)

# check how fantom annotates the mutations
data_icgc_wgs_to_join <- data_icgc_wgs %>%
  filter(!is.na(consequence_type)) %>% 
  select(chromosome = chromosome,
         start = chromosome_start,
         end = chromosome_end,
         icgc_mutation_id,
         icgc_donor_id,
         ref = reference_genome_allele,
         mut = mutated_to_allele) %>%
  distinct()

# number of mutation in each exon
data_exons_mutated <- data_exons_fantom %>% 
  genome_left_join(data_icgc_wgs_to_join) %>%
  select(chromosome = chromosome.x,
         start = start.x,
         end = end.x,
         icgc_mutation_id,
         icgc_donor_id) %>%
  mutate(count = ifelse(is.na(icgc_mutation_id), 0, 1)) %>%
  group_by(chromosome, start, end) %>%
  summarise(count = sum(count)) %>%
  mutate(length = end - start) %>%
  ungroup()

# exon mutations
mut_exon <- data_exons_fantom %>% 
  genome_left_join(data_icgc_wgs_to_join) %>%
  na.omit() %>%
  select(chromosome = chromosome.y,
         start = start.y,
         end = end.y,
         icgc_mutation_id,
         icgc_donor_id,
         ref,
         mut) %>%
  distinct()
mut_exon_bg <- mut_exon %>%
  mutate(start = start - 1,
         end = end + 1)

seq_exons_fantom <- getSeq(genome, names = data_exons_fantom$chromosome,
                               start = data_exons_fantom$start,
                               end = data_exons_fantom$end)
mat_tri <- reverseMerge(trinucleotideFrequency(seq_exons_fantom))
freq_tri <- enframe(colSums(mat_tri), name = "trinucleotide", value = "count") %>%
  mutate(percentage = 100*count/sum(count))

seq_exon_mutations <- getSeq(genome, names = mut_exon_bg$chromosome,
                                 start = mut_exon_bg$start,
                                 end = mut_exon_bg$end)
mat_tri_mut <- reverseMerge(trinucleotideFrequency(seq_exon_mutations))
freq_tri_mut <- enframe(colSums(mat_tri_mut), name = "trinucleotide", value = "mut_count") %>%
  arrange(desc(mut_count))

freq_tri <- freq_tri %>%
  inner_join(freq_tri_mut) %>%
  mutate(mut_rate = mut_count/count/num_donors)

data_exons_mutated <- data_exons_mutated %>%
  mutate(expected_count = as.vector(mat_tri %*% freq_tri$mut_rate * num_donors))

# annotate by refseq 
setwd(refseq.data.path)
mapping_hgnc <- read_delim("refseq_to_hgnc.txt", delim = "\t") %>%
  select(hgnc_symbol = `Approved symbol`,
         refseq_id = `RefSeq IDs`,
         refseq_id_ncbi = `RefSeq(supplied by NCBI)`) %>%
  pivot_longer(cols = refseq_id:refseq_id_ncbi, values_to = "gene") %>%
  select(hgnc_symbol, gene) %>%
  na.omit() %>%
  distinct()
mapping_refseq_exon <- read_delim("refseq_exons_171007.bed", delim = '\t',
                             col_names = F) %>%
  select(1:4) %>%
  rename(chromosome = "X1", start = "X2", end = "X3", gene = "X4") %>%
  mutate(gene = gsub("_exon.*$", "", gene))

data_exons_mutated_annotated <- data_exons_mutated %>%
  inner_join(mapping_refseq_exon) %>%
  inner_join(mapping_hgnc)

table_exon_mutation_by_gene <- data_exons_mutated_annotated %>% 
  group_by(hgnc_symbol) %>% 
  summarise(count = sum(count),
            exon_length = sum(length),
            expected_count = sum(expected_count)) %>% 
  arrange(desc(count)) 

