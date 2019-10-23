library(tidyverse)
library(fuzzyjoin)
library(knitr)
library(scales)
library(GenomicRanges)
options(tibble.width = Inf)

nonOverlap <- function(intervals){
  gr <- GRanges(seqnames = Rle(intervals$chromosome), 
                ranges = IRanges(intervals$start, intervals$end))
  gr_reduced <- GenomicRanges::reduce(gr)
  outcome <- tibble(
    chromosome = as.character(seqnames(gr_reduced)),
    start = start(gr_reduced),
    end = end(gr_reduced))
  return(outcome)
}

# file paths and file names
icgc.data.path <- "/Users/jz132/Desktop/Gordanlab/Data/ICGC"
genomic.interval.path <- "/Users/jz132/r_projects/cancer-mutations/pelinks"
filename <- "simple_somatic_mutation.open.LIAD-FR.tsv"

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

# import genomic coordinates of promoters, enhancers, and exons
setwd(genomic.interval.path)
data_enhancers_fantom <- read_delim("all_enhancers_fantom.bed", delim = "\t",
                                    col_names = c("chromosome", "start", "end"))
data_promoters_fantom <- read_delim("all_promoters_fantom.bed", delim = "\t",
                                    col_names = c("chromosome", "start", "end"))
data_exons_fantom <- read_delim("all_exons_fantom.bed", delim = "\t", 
                                col_names = c("chromosome", "start", "end"))
data_promoters_fantom <- nonOverlap(data_promoters_fantom)
data_exons_fantom <- nonOverlap(data_exons_fantom)
data_exons_in_promoters <- data_exons_fantom %>%
  genome_inner_join(data_promoters_fantom) %>%
  mutate(start = pmax(start.x, start.y),
         end = pmin(end.x, end.y)) %>%
  select(chromosome = chromosome.x, start, end) %>%
  distinct()

length_enhancers <- sum(data_enhancers_fantom$end - data_enhancers_fantom$start)
length_exons <- sum(data_exons_fantom$end - data_exons_fantom$start)
length_exons_in_promoters <- sum(data_exons_in_promoters$end - data_exons_in_promoters$start)
length_promoters <- sum(data_promoters_fantom$end - data_promoters_fantom$start) - length_exons_in_promoters

table_length <- tibble(position = c("promoter", "enhancer", "exon", "others"),
                       length = c(length_promoters, length_enhancers, length_exons, NA))

# check how fantom annotates the mutations
data_icgc_wgs_to_join <- data_icgc_wgs %>%
  filter(!is.na(consequence_type)) %>% 
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

mut_promoter <- data_promoters_fantom %>% 
  genome_left_join(data_icgc_wgs_to_join) %>%
  na.omit() %>%
  select(chromosome = chromosome.y,
         start = start.y,
         end = end.y,
         icgc_mutation_id,
         icgc_donor_id) %>%
  distinct() %>%
  dplyr::setdiff(mut_exon)

mut_others <- data_icgc_wgs_to_join %>%
  dplyr::setdiff(bind_rows(mut_exon, mut_promoter, mut_enhancer))

data_icgc_wgs_pos_fantom <- mutate(mut_exon, position = 'exon') %>%
  bind_rows(mutate(mut_promoter, position = 'promoter')) %>%
  bind_rows(mutate(mut_enhancer, position = 'enhancer')) %>%
  bind_rows(mutate(mut_others, position = 'others')) %>%
  select(icgc_mutation_id, icgc_donor_id, position)

table_mutation_pos_fantom <- data_icgc_wgs_pos_fantom %>% 
  group_by(position) %>%
  tally(name = "count") %>%
  inner_join(table_length) %>%
  mutate(mutation_rate = scientific(count/num_donors/length))

mut_enhancer_compare <- mut_enhancer %>%
  distinct(icgc_mutation_id) %>%
  inner_join(data_icgc_wgs %>% select(icgc_mutation_id, consequence_type)) %>%
  distinct()

mut_exon_compare <- mut_exon %>%
  distinct(icgc_mutation_id) %>%
  inner_join(data_icgc_wgs %>% select(icgc_mutation_id, consequence_type)) %>%
  distinct()

mut_promoter_compare <- mut_promoter %>%
  distinct(icgc_mutation_id) %>%
  inner_join(data_icgc_wgs %>% select(icgc_mutation_id, consequence_type)) %>%
  distinct()

mut_promoter_compare %>% 
  filter(consequence_type == "intron_variant" | 
           consequence_type == "upstream_gene_variant") %>%
  distinct(icgc_mutation_id)

mut_promoter_compare %>% 
  filter(consequence_type %in% nonexome) %>%
  distinct(icgc_mutation_id)

mut_enhancer_compare %>% 
  filter(consequence_type %in% nonexome) %>%
  distinct(icgc_mutation_id)
