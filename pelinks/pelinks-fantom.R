### In this code, we use the enhancer-tss association file
### from the Fantom project to build a promoter-enhancer links file. 
### This new file gives for each TSS, all the promoter and enhancer regions 
### associated with it.

# source packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("genomeIntervals")

library(tidyverse)
library(fuzzyjoin)
library(genomeIntervals)
library(GenomicRanges)
library(gtools)

rename <- dplyr::rename

# self-defined functions
isOverlap <- function(intervals, is.bed.format = F){
  gr <- GRanges(seqnames = Rle(intervals$chromosome), 
                ranges = IRanges(intervals$start + is.bed.format, intervals$end))
  
  print(ifelse(isDisjoint(gr), 
               "there is no overlap between intervals", 
               "there are overlaps between intervals"))
}

makeNonOverlap <- function(intervals, is.bed.format = F){
  gr <- GRanges(seqnames = Rle(intervals$chromosome), 
                ranges = IRanges(intervals$start + is.bed.format, intervals$end))
  gr_reduced <- GenomicRanges::reduce(gr)
  outcome <- tibble(
    chromosome = as.character(seqnames(gr_reduced)),
    start = start(gr_reduced),
    end = end(gr_reduced))
  return(outcome)
}

# input path
home.path <- "/Users/jz132/r_projects/cancer-mutations/pelinks" # set home directory
fantom.path <- paste0(home.path, "/FANTOM")
refseq.path <- paste0(home.path, "/RefSeq")

# output path
output.path <- home.path

# global parameters
options(stringsAsFactors = FALSE)
half_promoter <- 1000 # half promoter region length
half_tss <- half_promoter # enhancers within this distance to tss are removed
half_exon <- 200 # enhancers within this distance to exons are removed

## Part 1. Explore enhancer-promoter association and check the numbers in the paper
setwd(fantom.path)
data_enhancers_robust <- read_delim("robust_enhancers.bed", 
                                    delim = '\t', skip = 1, col_names = F) %>%
  select(chromosome = 1, start = 2, end = 3, enhancer = 4)
data_enhancers_permissive <- read_delim("permissive_enhancers.bed", 
                                        delim = '\t', skip = 1, col_names = F) %>%
  select(chromosome = 1, start = 2, end = 3, enhancer = 4)
data_association <-  read_delim("enhancer_tss_associations.bed", 
                                delim = '\t', skip = 1)

e_all_rob <- unique(data_enhancers_robust$enhancer) #all robust enhancers
e_all_per <- unique(data_enhancers_permissive$enhancer) #all permissive enhancers
all(e_all_rob %in% e_all_per) #robust enhancers are included in permissive ones

isOverlap(data_enhancers_permissive, is.bed.format = T) #check if enhancers overlap

refseq_info_split <- strsplit(as.character(data_association$name), split = ";")
e_eplinks_refseq <- unlist(lapply(refseq_info_split, `[[`, 1))
tss_eplinks_refseq <- unlist(lapply(refseq_info_split, `[[`, 2))
g_eplinks_refseq <- unlist(lapply(refseq_info_split, `[[`, 3))

all(e_eplinks_refseq %in% e_all_rob) #FALSE
all(e_eplinks_refseq %in% e_all_per) #TRUE 
# so eplinks_refseq is the enhancer-promoter links between 
# all permissive enhancers (not robust enhancers) and refseq TSS

# reproduce the numbers in the paper
length(tss_eplinks_refseq)/length(unique(e_eplinks_refseq)) #TSSs per enhancer
length(e_eplinks_refseq)/length(unique(tss_eplinks_refseq)) #enhancers per TSS
length(unique(e_eplinks_refseq))/length(e_all_per) #percentage of enhancers that have at least one TSS associated

## This is the end of using the original data ##
## It has several problems and is out of date ##


## Part 2. Filter out enhancers that are too close to exons or tss , and build tss-enhancer links
all_chromosomes <- unique(data_enhancers_permissive$chromosome)
setwd(refseq.path)
data_exons_refseq <- read_delim("refseq_exons_171007.bed", delim = "\t", col_names = F) %>%
  select(chromosome = 1, start = 2, end = 3) %>%
  filter(chromosome %in% all_chromosomes) %>%
  distinct()
data_exons_refseq_ext <- data_exons_refseq %>%
  select(chromosome, start, end) %>%
  mutate(start = start - half_exon, end = end + half_exon) %>%
  distinct()

data_tss_refseq <- read_delim("refseq_TSS_hg19_170929.bed", delim = '\t', col_names = F) %>%
  select(chromosome = 1, start = 2, end = 3, tss = 4, strand = 6) %>%
  filter(chromosome %in% all_chromosomes) %>%
  arrange(chromosome, start, end, tss) %>%
  distinct()
data_tss_refseq_ext <- data_tss_refseq %>%
  select(chromosome, start, end) %>%
  mutate(start = start - half_tss, end = end + half_tss) %>%
  distinct()

data_enhancers_exclude <- data_enhancers_permissive %>% 
  genome_inner_join(data_tss_refseq_ext) %>%
  select(chromosome.x, start.x, end.x) %>%
  bind_rows(data_enhancers_permissive %>% 
              genome_inner_join(data_exons_refseq_ext) %>%
              select(chromosome.x, start.x, end.x)) %>%
  distinct() %>%
  rename(chromosome = chromosome.x, start = start.x, end = end.x)

data_enhancers_new <- data_enhancers_permissive %>%
  anti_join(data_enhancers_exclude)
e_all_new <- unique(data_enhancers_new$enhancer)
list_drop <- which(!e_eplinks_refseq %in% e_all_new) # 4296 rows to exclude

# Update the enhancer, tss, and gene lists
refseq_info_split_filtered <- refseq_info_split[-list_drop]
e_eplinks_refseq <- unlist(lapply(refseq_info_split_filtered, `[[`, 1)) 
tss_eplinks_refseq <- unlist(lapply(refseq_info_split_filtered, `[[`, 2))
g_eplinks_refseq <- unlist(lapply(refseq_info_split_filtered, `[[`, 3))

# Get all the unique tss, note that some tss are no longer valid in RefSeq
tss_all <- unique(unlist(strsplit(tss_eplinks_refseq, split = ",")))
sum(!tss_all %in% data_tss_refseq$tss)  # 254 TSS no longer exist in the current RefSeq database

# Build tss-enhancer links (pelinks) from enhancer-tss links
# use valid tss in RefSeq only
eplinks_filtered <- tibble(enhancer = e_eplinks_refseq, 
                           tss = tss_eplinks_refseq, 
                           gene = g_eplinks_refseq) %>%
  mutate(tss = strsplit(tss, ",")) %>%
  unnest(tss) %>%
  filter(tss %in% data_tss_refseq$tss) %>%
  group_by(enhancer, gene) %>%
  summarise(tss = paste0(tss, collapse=";")) %>%
  select(enhancer, tss, gene) %>%
  ungroup()

mapping_tss_enhancer <- eplinks_filtered %>% 
  select(tss, enhancer) %>%
  group_by(tss) %>%
  summarise(enhancer = paste0(enhancer, collapse=";")) %>% 
  ungroup()

mapping_tss_gene <- eplinks_filtered %>% 
  select(tss, gene) %>%
  distinct()

# For each TSS, find all enhancers linked to it, and denote by NA if none
pelinks_raw_unnest <- mapping_tss_enhancer %>%
  mutate(tss = strsplit(tss, ";")) %>%
  unnest(tss)

data_tss_refseq <- data_tss_refseq %>%
  select(chromosome, tss_pos = start, tss, strand)

pelinks_sameTSS <- data_tss_refseq %>% 
  left_join(pelinks_raw_unnest, by = "tss") %>%
  mutate(e_count = ifelse(is.na(enhancer), 0, sapply(strsplit(enhancer, ";"), length)))


## Part 3. Define raw promoters as the regions close to TSS and add it to pelinks
data_promoters_raw <- pelinks_sameTSS %>% 
  select(chromosome, tss_pos, tss, strand) %>% 
  mutate(start = tss_pos - half_promoter*ifelse(strand == "+", 1, -1),
         end = tss_pos + half_promoter*ifelse(strand == "+", 1, -1)) %>%
  select(-tss_pos)

# Add the promoter variable in pelinks_sameTSS
pelinks_sameTSS_complete <- pelinks_sameTSS %>% 
  bind_cols(data_promoters_raw %>% select(start, end)) %>%
  mutate(promoter = paste0(chromosome, ":", start, "-", end)) %>%
  group_by(chromosome, tss, tss_pos, enhancer, e_count) %>%
  summarise(promoter = paste0(promoter, collapse = ";")) %>%
  ungroup() %>%
  select(tss, chromosome, tss_pos, promoter, enhancer, e_count)

## Output the result: promoter-enhancer links, all promoters, all enhancers, all exons
data_enhancers_output <- data_enhancers_new
data_promoters_output <- data_promoters_raw %>% 
  select(chromosome, start, end, strand) %>%
  mutate(promoter = paste0(chromosome, ":", start, "-", end)) %>% 
  select(chromosome, start, end, promoter, strand)
data_exons_output <- data_exons_refseq %>% 
  mutate(exon = paste0(chromosome, ":", start, "-", end))

setwd(output.path)
write.csv(eplinks_filtered, "eplinks-fantom-filtered.csv", quote = F, row.names = F)
write.csv(pelinks_sameTSS_complete, "pelinks-fantom-filtered.csv", quote = F, row.names = F)
write.table(data_enhancers_output, "all_enhancers_fantom.txt", 
            quote = F, row.names = F, col.names = F, sep = "\t")
write.table(data_promoters_output, "all_promoters_refseq.txt", 
            quote = F, row.names = F, col.names = F, sep = "\t")
write.table(data_exons_output, "all_exons_refseq.txt", 
            quote = F, row.names = F, col.names = F, sep = "\t")

# # a histogram showing the distribution of enhancer length
# setwd(output.path)
# data_enhancers_plot <- data_enhancers_permissive %>%
#   mutate(length = end - start)
# png(file = "enhancer_length_dist.png", width = 1200, height = 800, res = 160)
# ggplot(data_enhancers_plot, aes(x = length)) +
#   geom_histogram(fill = "lightblue", alpha = 0.7, boundary = 0) +
#   labs(x = "enhancer length") +
#   theme_bw()
# dev.off()
# 
# # A bar plot showing the number of enhancers each tss is associated with
# setwd(output.path)
# pelinks_plot <- pelinks_sameTSS %>%
#   mutate(category = cut(e_count, breaks = c(0, 5, 10, 15, 20, Inf), right = F)) %>%
#   group_by(category) %>%
#   summarise(count = n())
# png(file = "enhancers_per_tss.png", width = 1200, height = 800, res = 160)
# ggplot(pelinks_plot, aes(x = category, y = count)) +
#   geom_bar(stat = "identity", fill = "lightblue", alpha = 0.7, width = 0.5) + 
#   geom_text(aes(label = count), nudge_y = 500) +
#   labs(x = "number of enhancers per tss") +
#   theme_bw()
# dev.off()

