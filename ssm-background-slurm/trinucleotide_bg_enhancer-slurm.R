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
icgc.data.path <- "/data/gordanlab/jingkang/Desktop/Gordanlab/Data/ICGC"
genomic.interval.path <- "~/r_projects/cancer-mutations/pelinks"
output.path <- "~/r_projects/cancer-mutations/ssm-background-slurm/"

datasetname <- "LIRI-JP"
filename <- paste0("simple_somatic_mutation.open.", datasetname, ".tsv")

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

# focus on single base substitution in WGS only
data_icgc_wgs <- data_icgc_raw %>% 
  filter(sequencing_strategy == "WGS", mutation_type == "single base substitution")
num_donors <- length(data_icgc_wgs %>% distinct(icgc_donor_id) %>% pull())
num_mutations <- length(data_icgc_wgs %>% distinct(icgc_mutation_id, icgc_donor_id) %>% pull())
if(!all(data_icgc_wgs$mutated_from_allele == data_icgc_wgs$reference_genome_allele)){
  warning("reference genome allele is not the same as mutated from allele")
}

# import genomic coordinates of enhancers and links between enhancers and tss
setwd(genomic.interval.path)
data_enhancers_fantom <- read_delim("all_enhancers_fantom.txt", delim = "\t",
                                    col_names = c("chromosome", "start", "end", "enhancer")) %>%
  arrange(chromosome, start, end)

data_exons_refseq <- read_delim("all_exons_refseq.txt", delim = "\t", 
                                col_names = c("chromosome", "start", "end", "exon"))

# count the number of mutations in each promoter
data_icgc_wgs_to_join <- data_icgc_wgs %>%
  filter(!is.na(consequence_type)) %>% 
  select(chromosome = chromosome,
         start = chromosome_start,
         end = chromosome_end,
         icgc_mutation_id,
         icgc_donor_id,
         ref = reference_genome_allele,
         mut = mutated_to_allele) %>%
  genome_anti_join(data_exons_refseq, by = c("chromosome", "start", "end")) %>%
  distinct()

data_enhancers_mutated <- data_enhancers_fantom %>% 
  genome_left_join(data_icgc_wgs_to_join) %>%
  select(chromosome = chromosome.x,
         start = start.x,
         end = end.x,
         enhancer,
         icgc_mutation_id,
         icgc_donor_id) %>%
  mutate(count = ifelse(is.na(icgc_mutation_id), 0, 1)) %>%
  group_by(chromosome, start, end, enhancer) %>%
  summarise(count = sum(count)) %>%
  mutate(length = end - start + 1) %>%
  ungroup()

# enhancer mutations
mut_enhancer <- data_enhancers_fantom %>% 
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

# enhancer sequences and trinucleotide frequencies
seq_enhancers_fantom <- getSeq(genome, names = data_enhancers_fantom$chromosome,
                               start = data_enhancers_fantom$start,
                               end = data_enhancers_fantom$end)
mat_tri <- reverseMerge(trinucleotideFrequency(seq_enhancers_fantom))
freq_tri <- enframe(colSums(mat_tri), name = "trinucleotide", value = "count")

# trinucleotide background of the mutations
mut_enhancer_bg <- mut_enhancer %>%
  mutate(start = start - 1,
         end = end + 1)

# frequencies of trinucleotides that are mutated
seq_enhancer_mutations_ref <- getSeq(genome, names = mut_enhancer_bg$chromosome,
                                     start = mut_enhancer_bg$start,
                                     end = mut_enhancer_bg$end)
seq_enhancer_mutations_mut <- replaceLetterAt(seq_enhancer_mutations_ref, 
                                              at = matrix(c(F, T, F), 
                                                          nrow = length(seq_enhancer_mutations_ref),
                                                          ncol = 3,
                                                          byrow = T), 
                                              mut_enhancer_bg$mut)

table_mutation_tri <- tibble(
  ref = as.character(seq_enhancer_mutations_ref),
  ref_rev = as.character(reverseComplement(seq_enhancer_mutations_ref)),
  mut = as.character(seq_enhancer_mutations_mut),
  mut_rev =  as.character(reverseComplement(seq_enhancer_mutations_mut)),
) %>%
  mutate(ref = ifelse(ref < ref_rev, ref, ref_rev),
         mut = ifelse(ref < ref_rev, mut, mut_rev)) %>%
  group_by(ref, mut) %>%
  tally(name = "count") %>% 
  ungroup()

mat_tri_mut <- reverseMerge(trinucleotideFrequency(seq_enhancer_mutations_ref))
freq_tri_mut <- enframe(colSums(mat_tri_mut), name = "trinucleotide", value = "mut_count") %>%
  arrange(desc(mut_count))

freq_tri <- freq_tri %>%
  inner_join(freq_tri_mut) %>%
  mutate(mut_rate = mut_count/count/num_donors)

print("data_enhancers_mutated created")


