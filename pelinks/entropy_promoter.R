library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
library(fuzzyjoin)

options(tibble.width = Inf)
slice <- dplyr::slice
rename <- dplyr::rename

### start functions ###

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

## function DNAToBin ##
## transform DNA to binary number
## to further transform binary to decimal number, use base R function strtoi
DNAToBin <- function(DNA){
  temp <- DNA
  temp <- gsub("A", "00", temp)
  temp <- gsub("C", "01", temp)
  temp <- gsub("G", "10", temp)
  temp <- gsub("T", "11", temp)
  return(temp)
}

## function absmax ##
## get the element with maximum absolute value from a vector of numbers
absmax <- function(x) {x[which.max(abs(x))]}

### end functions ###

# file paths and file names
icgc.data.path <- "~/Desktop/Gordanlab/Data/ICGC"
genomic.interval.path <- "~/r_projects/cancer-mutations/pelinks"
output.path <- "~/r_projects/cancer-mutations/ssm-background/output"

datasetname <- "PACA-CA"
filename_icgc_ssm <- paste0("simple_somatic_mutation.open.", datasetname, ".tsv")
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
data_icgc_try <- read_tsv(filename_icgc_ssm, col_names = T, n_max = 5)
col_indicator <- paste(ifelse(colnames(data_icgc_try) %in% keep_cols, "?", "-"), collapse = "")
data_icgc_raw <- read_tsv(filename_icgc_ssm, col_names = T, col_types = col_indicator) %>%
  mutate(chromosome = paste0("chr", chromosome))

# focus on single base substitution in WGS only
data_icgc_wgs <- data_icgc_raw %>% 
  filter(sequencing_strategy == "WGS", mutation_type == "single base substitution")
num_donors <- length(data_icgc_wgs %>% distinct(icgc_donor_id) %>% pull())
num_mutations <- length(data_icgc_wgs %>% distinct(icgc_mutation_id, icgc_donor_id) %>% pull())
if(!all(data_icgc_wgs$mutated_from_allele == data_icgc_wgs$reference_genome_allele)){
  warning("reference genome allele is not the same as mutated from allele")
}

# import genomic coordinates of promoters
setwd(genomic.interval.path)
data_promoters_refseq <- read_delim("all_promoters_refseq.txt", delim = "\t",
                                    col_names = c("chromosome", "start", "end", "promoter")) %>%
  arrange(chromosome, start, end) %>% 
  distinct()

data_exons_refseq <- read_delim("all_exons_refseq.txt", delim = "\t", 
                                col_names = c("chromosome", "start", "end", "exon"))

# keep only promoters for protein coding genes
mapping_refseq_tss <- read_delim("refseq_TSS_hg19_170929.bed", delim = '\t',
                                 col_names = F) %>%
  select(c(1,2,4)) %>%
  dplyr::rename(chromosome = "X1", pos = "X2", tss = "X4") %>% 
  filter(grepl("NM", tss))

mapping_hgnc <- read_delim("refseq_to_hgnc.txt", delim = "\t") %>%
  select(gene = `Approved symbol`,
         refseq_id = `RefSeq IDs`,
         refseq_id_ncbi = `RefSeq(supplied by NCBI)`) %>%
  pivot_longer(cols = refseq_id:refseq_id_ncbi, values_to = "tss") %>%
  select(gene, tss) %>%
  filter(grepl("NM", tss)) %>%
  na.omit() %>%
  distinct()

data_promoters_refseq <- data_promoters_refseq %>% 
  mutate(pos = (start + end)/2) %>% 
  inner_join(mapping_refseq_tss) %>%
  inner_join(mapping_hgnc) %>%
  select(-c("pos", "tss", "gene")) %>% 
  distinct()

# filter out mutations that are not in promoters
data_mutations_promoter <- data_icgc_wgs %>%
  filter(!is.na(consequence_type)) %>% 
  select(chromosome = chromosome,
         start = chromosome_start,
         end = chromosome_end,
         icgc_mutation_id,
         icgc_donor_id,
         ref = reference_genome_allele,
         mut = mutated_to_allele) %>%
  genome_anti_join(data_exons_refseq, by = c("chromosome", "start", "end")) %>%
  genome_inner_join(data_promoters_refseq, by = c("chromosome", "start", "end")) %>%
  select(chromosome = chromosome.x,
         start = start.x,
         end = end.x,
         icgc_mutation_id, icgc_donor_id, ref, mut) %>%
  distinct()

# promoter sequences and trinucleotide frequencies
seq_promoters_refseq <- getSeq(genome, names = data_promoters_refseq$chromosome,
                               start = data_promoters_refseq$start,
                               end = data_promoters_refseq$end)
seq_promoters_refseq_bg <- getSeq(genome, names = data_promoters_refseq$chromosome,
                                  start = data_promoters_refseq$start - 1,
                                  end = data_promoters_refseq$end + 1)
mat_tri <- reverseMerge(trinucleotideFrequency(seq_promoters_refseq_bg))
freq_tri <- enframe(colSums(mat_tri), name = "trinucleotide", value = "count")

# frequencies of trinucleotides that are mutated
seq_promoter_mutations_ref <- getSeq(genome, names = data_mutations_promoter$chromosome,
                                     start = data_mutations_promoter$start - 1,
                                     end = data_mutations_promoter$end + 1)
seq_promoter_mutations_mut <- replaceLetterAt(seq_promoter_mutations_ref, 
                                              at = matrix(c(F, T, F), 
                                                          nrow = length(seq_promoter_mutations_ref),
                                                          ncol = 3,
                                                          byrow = T), 
                                              data_mutations_promoter$mut)

table_mutation_tri <- tibble(
  ref = as.character(seq_promoter_mutations_ref),
  ref_rev = as.character(reverseComplement(seq_promoter_mutations_ref)),
  mut = as.character(seq_promoter_mutations_mut),
  mut_rev =  as.character(reverseComplement(seq_promoter_mutations_mut)),
) %>%
  mutate(ref = ifelse(ref < ref_rev, ref, ref_rev),
         mut = ifelse(ref < ref_rev, mut, mut_rev)) %>%
  group_by(ref, mut) %>%
  tally(name = "count") %>% 
  ungroup()

mat_tri_mut <- reverseMerge(trinucleotideFrequency(seq_promoter_mutations_ref))
freq_tri_mut <- enframe(colSums(mat_tri_mut), name = "trinucleotide", value = "mut_count")

freq_tri_mut_rate <- freq_tri_mut %>%
  inner_join(freq_tri, by = "trinucleotide") %>%
  mutate(mut_rate = mut_count/count/num_donors)

table_mutation_tri_mut_rate <- table_mutation_tri %>% 
  group_by(ref) %>% 
  mutate(prop = count/sum(count)) %>% 
  ungroup() %>% 
  inner_join(freq_tri_mut_rate %>% select(ref = trinucleotide, mut_rate), by = "ref") %>%
  mutate(tri_mut_rate = prop*mut_rate) %>% 
  select(ref_tri = ref,
         mut_tri = mut,
         tri_mut_rate) %>% 
  mutate(mut_type = row_number())

# for each promoter, count the number of mutations in it
data_promoters_mutated <- data_promoters_refseq %>% 
  genome_left_join(data_mutations_promoter, by = c("chromosome", "start", "end")) %>%
  select(chromosome = chromosome.x,
         start = start.x,
         end = end.x,
         promoter,
         icgc_mutation_id,
         icgc_donor_id) %>%
  mutate(count = ifelse(is.na(icgc_mutation_id), 0, 1)) %>%
  group_by(chromosome, start, end, promoter) %>%
  summarise(count = sum(count)) %>%
  mutate(length = end - start + 1) %>%
  ungroup()

promoters_w_mutation <- which(data_promoters_mutated$count > 0)
n_promoters_w_mutation <- length(promoters_w_mutation)

# calculate Shannon entropy for each mutated promoter

s_entropy <- NULL

for(i in 1:n_promoters_w_mutation){
  promoter_ex <- promoters_w_mutation[i]
  print(promoter_ex)
  
  data_promoter_ex <- data_promoters_refseq %>% slice(promoter_ex) # the promoter location
  seq_promoter_ex <- seq_promoters_refseq[[promoter_ex]] # the promoter sequence
  seq_promoter_ex_bg <- getSeq(genome, data_promoter_ex$chromosome, 
                               data_promoter_ex$start - 5, 
                               data_promoter_ex$end + 5) # the promoter sequence including 5bp context
  seq_promoter_ex_11mer <- DNAStringSet(seq_promoter_ex_bg, 
                                        start = 1:(length(seq_promoter_ex_bg)-10),
                                        width = 11) # cut the promoter into 11mers
  
  # generate all possible single mutations
  data_promoter_all_possible_mutations <- tibble(
    chromosome = data_promoter_ex$chromosome,
    pos = seq(data_promoter_ex$start, data_promoter_ex$end),
    promoter = promoter_ex,
    seq_mut_bg = as.character(seq_promoter_ex_11mer)
  ) %>% 
    mutate(ref = unlist(strsplit(as.character(seq_promoter_ex), ""))) %>% 
    slice(rep(seq(1, length(seq_promoter_ex)), each = 4)) %>% 
    mutate(mut = rep(c("A", "C", "G", "T"), length(seq_promoter_ex))) %>% 
    filter(ref != mut) %>% 
    mutate(twimer = paste0(seq_mut_bg, mut)) %>%
    mutate(idx = DNAToBin(twimer)) %>%
    mutate(idx = strtoi(idx, base = 2)) %>%
    select(-seq_mut_bg) %>% 
    mutate(ref_tri = substr(twimer, 5, 7)) %>%
    mutate(mut_tri = paste0(substr(ref_tri, 1, 1), mut, substr(ref_tri, 3, 3))) %>% 
    mutate(ref_tri_rev = as.character(reverseComplement(DNAStringSet(ref_tri))),
           mut_tri_rev = as.character(reverseComplement(DNAStringSet(mut_tri)))) %>% 
    mutate(ref_tri = ifelse(ref_tri < ref_tri_rev, ref_tri, ref_tri_rev),
           mut_tri = ifelse(ref_tri < ref_tri_rev, mut_tri, mut_tri_rev)) %>% 
    select(-c(ref_tri_rev, mut_tri_rev))  %>% 
    inner_join(table_mutation_tri_mut_rate, by = c("ref_tri", "mut_tri")) %>% 
    mutate(cond_tri_mut_rate = tri_mut_rate/sum(tri_mut_rate))
  
  data_promoter_mutations_by_twimer <- data_promoter_all_possible_mutations %>%
    group_by(twimer, ref_tri, mut_tri) %>% 
    summarise(prob = sum(cond_tri_mut_rate)) %>% 
    ungroup()
  
  s_entropy_new <- sum(-data_promoter_mutations_by_twimer$prob*
                         log2(data_promoter_mutations_by_twimer$prob))
  s_entropy <- c(s_entropy, s_entropy_new)
  
  print(s_entropy_new)
}

data_promoters_entropy <- data_promoters_mutated %>%
  filter(count > 0) %>% 
  mutate(s_entropy = s_entropy)

setwd(output.path)
filename_output <- paste0("promoters_entropy_", datasetname, ".csv")
write_csv(data_promoters_entropy, filename_output)

