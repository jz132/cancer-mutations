library(tidyverse)
library(fuzzyjoin)
options(tibble.width = Inf)

# file path
icgc.data.path <- "/Users/jz132/Desktop/Gordanlab/Data/ICGC"
genomic.interval.path <- "/Users/jz132/r_projects/cancer-mutations/pelinks"
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

# import ICGC cancer mutations
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

# focus on single base substitution in WGS only
data_icgc_wgs <- data_icgc_raw %>% 
  filter(sequencing_strategy == "WGS", mutation_type == "single base substitution")

rm(data_icgc_raw)
data_icgc_wgs_analysis <- data_icgc_wgs %>%
  distinct(icgc_mutation_id, icgc_donor_id, 
           chromosome, chromosome_start, chromosome_end,
           reference_genome_allele, mutated_from_allele, mutated_to_allele,
           consequence_type) %>%
  filter(!is.na(consequence_type)) %>% 
  mutate(position = ifelse(consequence_type %in% exome, "exome", "nonexome"))

table_mutation_consequence_snpeff <- data_icgc_wgs_analysis %>%
  group_by(consequence_type) %>%
  tally()

# check how snpeff annotates the mutations
data_icgc_wgs_pos_snpeff <- data_icgc_wgs_analysis %>%
  group_by(icgc_mutation_id, icgc_donor_id) %>%
  summarise(position = case_when(
    all(position == "exome") ~ "exome",
    all(position == "nonexome") ~ "nonexome",
    TRUE ~ "both"
  )) %>%
  ungroup()

table_mutation_pos_snpeff <- data_icgc_wgs_pos_snpeff %>% 
  group_by(position) %>%
  tally()

# import genomic coordinates of promoters, enhancers, and exons
setwd(genomic.interval.path)
data_promoters_fantom <- read_delim("all_promoters_fantom.bed", delim = "\t",
                                    col_names = c("chromosome", "start", "end"))
data_enhancers_fantom <- read_delim("all_enhancers_fantom.bed", delim = "\t",
                                    col_names = c("chromosome", "start", "end"))
data_exons_fantom <- read_delim("all_exons_fantom.bed", delim = "\t", 
                                col_names = c("chromosome", "start", "end"))

# check how fantom annotates the mutations
data_icgc_wgs_to_join <- data_icgc_wgs_analysis %>%
  select(chromosome = chromosome,
         start = chromosome_start,
         end = chromosome_end,
         icgc_mutation_id,
         icgc_donor_id) %>%
  distinct()

mut_enhancer <- data_enhancers_fantom %>% 
  genome_left_join(data_icgc_wgs_to_join) %>%
  na.omit() %>%
  select(chromosome = chromosome.y,
         start = start.y,
         end = end.y,
         icgc_mutation_id,
         icgc_donor_id) %>%
  distinct()

mut_exon <- data_exons_fantom %>% 
  genome_left_join(data_icgc_wgs_to_join) %>%
  na.omit() %>%
  select(chromosome = chromosome.y,
         start = start.y,
         end = end.y,
         icgc_mutation_id,
         icgc_donor_id) %>%
  distinct()

mut_promoter_raw <- data_promoters_fantom %>% 
  genome_left_join(data_icgc_wgs_to_join) %>%
  na.omit() %>%
  select(chromosome = chromosome.y,
         start = start.y,
         end = end.y,
         icgc_mutation_id,
         icgc_donor_id) %>%
  distinct()

mut_promoter <- mut_promoter_raw %>%
  setdiff(mut_exon)

mut_others <- data_icgc_wgs_to_join %>%
  setdiff(bind_rows(mut_exon, mut_promoter, mut_enhancer))

data_icgc_wgs_pos_fantom <- mutate(mut_exon, position = 'exon') %>%
  bind_rows(mutate(mut_promoter, position = 'promoter')) %>%
  bind_rows(mutate(mut_enhancer, position = 'enhancer')) %>%
  bind_rows(mutate(mut_others, position = 'others')) %>%
  select(icgc_mutation_id, icgc_donor_id, position)

table_mutation_pos_fantom <- data_icgc_wgs_pos_fantom %>% 
  group_by(position) %>%
  tally()

# output the result
data_icgc_wgs_output <- data_icgc_wgs %>%
  select(chromosome = chromosome,
         start = chromosome_start,
         end = chromosome_end,
         ref = mutated_from_allele,
         alt = mutated_to_allele)

setwd(output.path)
write.table(data_icgc_wgs_output, "single_mutations_LIRI-JP.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")


