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
figure.path <- "/Users/jz132/r_projects/cancer-mutations/ssm-background/Figures/"
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

# import genomic coordinates of promoters
setwd(genomic.interval.path)
data_promoters_refseq <- read_delim("all_promoters_refseq.txt", delim = "\t", 
                                col_names = c("chromosome", "start", "end", "promoter")) %>%
  arrange(chromosome, start, end) %>%
  distinct()

data_promoters_refseq <- data_promoters_refseq %>%
  mutate(end = end - 1000) # for now, we only consider the part upstream of tss

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
  distinct()

data_promoters_mutated <- data_promoters_refseq %>% 
  genome_left_join(data_icgc_wgs_to_join) %>%
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

# promoter mutations
mut_promoter <- data_promoters_refseq %>% 
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

# promoter sequences and trinucleotide frequencies
seq_promoters_refseq <- getSeq(genome, names = data_promoters_refseq$chromosome,
                               start = data_promoters_refseq$start,
                               end = data_promoters_refseq$end)
mat_tri <- reverseMerge(trinucleotideFrequency(seq_promoters_refseq))
freq_tri <- enframe(colSums(mat_tri), name = "trinucleotide", value = "count")

setwd(figure.path)
png(file = "promoter_trinucleotide_dist.png", width = 1200, height = 800, res = 160)
ggplot(data = freq_tri, aes(x = trinucleotide, y = count)) + 
  geom_col() + 
  ggtitle("LIRI-JP promoter sequence") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45))
dev.off()

# trinucleotide background of the mutations
mut_promoter_bg <- mut_promoter %>%
  mutate(start = start - 1,
         end = end + 1)

# frequencies of trinucleotides that are mutated
seq_promoter_mutations_ref <- getSeq(genome, names = mut_promoter_bg$chromosome,
                                     start = mut_promoter_bg$start,
                                     end = mut_promoter_bg$end)
seq_promoter_mutations_mut <- replaceLetterAt(seq_promoter_mutations_ref, 
                                              at = matrix(c(F, T, F), 
                                                          nrow = length(seq_promoter_mutations_ref),
                                                          ncol = 3,
                                                          byrow = T), 
                                              mut_promoter_bg$mut)
table_mutation_tri <- tibble(
  ref = as.character(seq_promoter_mutations_ref),
  ref_rev = as.character(reverseComplement(seq_promoter_mutations_ref)),
  mut = as.character(seq_promoter_mutations_mut),
  mut_rev =  as.character(reverseComplement(seq_promoter_mutations_mut)),
) %>%
  mutate(ref = ifelse(ref < ref_rev, ref, ref_rev),
         mut = ifelse(ref < ref_rev, mut, mut_rev)) %>%
  group_by(ref, mut) %>%
  tally(name = "count")

setwd(figure.path)
png(file = "promoter_mutated_trinucleotide_dist.png", width = 1200, height = 800, res = 160)
ggplot(data = table_mutation_tri, mapping = aes(x = ref, y = count, 
                                                fill = paste0(substr(ref, 2, 2),
                                                              "->",
                                                              substr(mut, 2, 2)))) + 
  geom_col() + 
  labs(fill = "Mutation") +
  ggtitle("LIRI-JP promoter mutation") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45))
dev.off()

mat_tri_mut <- reverseMerge(trinucleotideFrequency(seq_promoter_mutations_ref))
freq_tri_mut <- enframe(colSums(mat_tri_mut), name = "trinucleotide", value = "mut_count") %>%
  arrange(desc(mut_count))

freq_tri <- freq_tri %>%
  inner_join(freq_tri_mut) %>%
  mutate(mut_rate = mut_count/count/num_donors)

data_promoters_mutated <- data_promoters_mutated %>%
  mutate(expected_count = as.vector(mat_tri %*% freq_tri$mut_rate * num_donors),
         var_count = as.vector(mat_tri %*% 
                                 freq_tri$mut_rate*(1-freq_tri$mut_rate) * num_donors))

# annotate by refseq 
setwd(refseq.data.path)
mapping_hgnc <- read_delim("refseq_to_hgnc.txt", delim = "\t") %>%
  select(gene = `Approved symbol`,
         refseq_id = `RefSeq IDs`,
         refseq_id_ncbi = `RefSeq(supplied by NCBI)`) %>%
  pivot_longer(cols = refseq_id:refseq_id_ncbi, values_to = "tss") %>%
  select(gene, tss) %>%
  na.omit() %>%
  distinct()
mapping_refseq_tss <- read_delim("refseq_TSS_hg19_170929.bed", delim = '\t',
                             col_names = F) %>%
  select(c(1,2,4)) %>%
  dplyr::rename(chromosome = "X1", pos = "X2", tss = "X4")

data_promoters_mutated_annotated <- data_promoters_mutated %>%
  mutate(pos = end) %>%
  inner_join(mapping_refseq_tss) %>%
  inner_join(mapping_hgnc) %>%
  select(-pos)

# per gene
table_promoter_mutation_by_gene <- data_promoters_mutated_annotated %>%
  group_by(gene) %>%
  summarise(count = sum(count),
            promoter_length = sum(length),
            expected_count = sum(expected_count),
            var_count = sum(var_count)) %>%
  mutate(p_value = ppois(count, expected_count, lower.tail = F)) %>%
  arrange(desc(count))

# per tss
table_promoter_mutation_by_tss <- data_promoters_mutated_annotated %>%
  group_by(tss) %>%
  summarise(count = sum(count),
            promoter_length = sum(length),
            expected_count = sum(expected_count),
            var_count = sum(var_count)) %>%
  mutate(p_value = ppois(count, expected_count, lower.tail = F)) %>%
  arrange(desc(count))

