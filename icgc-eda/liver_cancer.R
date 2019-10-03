library(tidyverse)
options(tibble.width = Inf)

# file path
icgc.data.path <- "/Users/jz132/Desktop/Gordanlab/Data/ICGC/"
output.path <- "/Users/jz132/Desktop/Gordanlab/Main_Quests/icgc_data_analysis/output"

# global parameters
appr_genome_len <- 3*10^9
appr_exome_len <- 3*10^7

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

# import ICGC breast cancer mutations
# data source https://dcc.icgc.org/releases/current/Projects/
setwd(icgc.data.path)
filename <- "simple_somatic_mutation.open.LIRI-JP.tsv"
data_icgc_try <- read_tsv(filename, col_names = T, n_max = 5)
col_indicator <- paste(ifelse(colnames(data_icgc_try) %in% keep_cols, "?", "-"), 
                       collapse = "")

data_icgc_raw <- read_tsv(filename,
                          col_names = T, col_types = col_indicator) %>%
  mutate(chromosome = paste0("chr", chromosome))

all_types <- data_icgc_raw %>% distinct(consequence_type) %>% pull()
nonexome <- c("intron_variant", "intragenic_variant", "upstream_gene_variant", 
              "downstream_gene_variant", "intergenic_region")
exome <- setdiff(all_types, c(nonexome, NA))

# reproduce some numbers on the icgc website
data_icgc_raw %>% distinct(icgc_donor_id) %>% nrow() # total number of donors
data_icgc_raw %>% 
  distinct(icgc_donor_id ,icgc_mutation_id) %>% nrow() # total number of mutations
data_icgc_raw %>% filter(consequence_type %in% exome) %>% 
  distinct(icgc_donor_id, icgc_mutation_id) %>% nrow() # total number of exome mutations

# median number of exome mutations per mb
data_icgc_in_exome_by_donor <- data_icgc_raw %>%
  filter(consequence_type %in% exome) %>%
  group_by(icgc_donor_id) %>%
  summarise(n_mutation = n_distinct(icgc_mutation_id))
median(data_icgc_in_exome_by_donor$n_mutation)*10^6/appr_exome_len 

# focus on single base substitution and analyze wxs and wgs separately
data_icgc_wxs <- data_icgc_raw %>% 
  filter(sequencing_strategy == "WXS", mutation_type == "single base substitution")
data_icgc_wgs <- data_icgc_raw %>% 
  filter(sequencing_strategy == "WGS", mutation_type == "single base substitution")

rm(data_icgc_raw)
data_icgc_wgs_analysis <- data_icgc_wgs %>%
  distinct(icgc_mutation_id, icgc_donor_id, 
           chromosome, chromosome_start, chromosome_end,
           reference_genome_allele, mutated_from_allele, mutated_to_allele,
           consequence_type) %>%
  mutate(position = ifelse(consequence_type %in% exome, "exome", "nonexome"))

data_icgc_wgs_position <- data_icgc_wgs_analysis %>%
  group_by(icgc_mutation_id, icgc_donor_id) %>%
  summarise(position = case_when(
    all(position == "exome") ~ "exome",
    all(position == "nonexome") ~ "nonexome",
    TRUE ~ "both"
  )) %>%
  ungroup()

data_icgc_wgs_position %>% 
  group_by(position) %>%
  tally()
