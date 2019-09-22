### In this code, we use primarily the enhancer-tss association file
### from the Fantom project to build a pelinks_sameTSS file. 
### This new file gives for each TSS, all the promoter and enhancer regions 
### associatited with it.

# source packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("genomeIntervals")
library(genomeIntervals)
library(tidyverse)
library(dplyr)
library(tidyr)
library(fuzzyjoin)
library(gtools)

# input path
home.path <- getwd()
fantom.path <- paste0(home.path, "/FANTOM")
refseq.path <- paste0(home.path, "/RefSeq")
processed.path <- paste0(home.path, "/Processed")

# output path
output.path <- home.path

# global parameters
options(stringsAsFactors = FALSE)
half_p <- 2000 # half promoter region length
half_tss <- 500 # enhancers within this distance to tss are removed
half_exon <- 200 # enhancers within this distance to exons are removed

## Part 1. Explore enhancer-promoter association and check the numbers in the paper
setwd(fantom.path)
data_enhancers_rob <- as_tibble(read.table("robust_enhancers.bed", 
                                           header = F, sep = '\t')) %>%
  select(V1, V2, V3, V4) %>%
  rename(chromosome = V1, start = V2, end = V3, enhancer = V4)
data_enhancers_per <- as_tibble(read.table("permissive_enhancers.bed", 
                                           header = F, sep = '\t', skip = 1)) %>%
  select(V1, V2, V3, V4) %>%
  rename(chromosome = V1, start = V2, end = V3, enhancer = V4)
data_association_refseq <- as_tibble(read.table("enhancer_tss_associations_complete.bed", 
                                      header = F, sep = '\t'))

e_all_rob <- unique(data_enhancers_rob$enhancer) #all robust enhancers
e_all_per <- unique(data_enhancers_per$enhancer) #all permissive enhancers
all(e_all_rob %in% e_all_per) #robust enhancers are included in permissive ones

# a histogram showing the distribution of enhancer length
setwd(output.path)
data_enhancers_plot <- data_enhancers_per %>%
  mutate(length = end - start)
png(file = "e_length.png", width = 1200, height = 800, res = 160)
ggplot(data_enhancers_plot, aes(x = length)) +
  geom_histogram(fill = "lightblue", alpha = 0.7, boundary = 0) +
  labs(x = "enhancer length") +
  theme_bw()
dev.off()
# about 0.4% of the whole genome that are enhancers defined by the paper
sum(data_enhancers_plot %>% pull(length))/(3*10^9) 

enhancers_check_result <- data_enhancers_per %>% 
  genome_inner_join(data_enhancers_per, 
                    by = c("chromosome", "start", "end")) %>% 
  filter(start.x != start.y | end.x != end.y) %>% #exclude self
  filter(start.x != end.y & start.y != end.x) #due to bed file interval format (,]
print(ifelse(nrow(enhancers_check_result) == 0, 
             "no overlap between enhancers", "there is overlap between enhancers"))

refseq_info_split <- strsplit(as.character(data_association_refseq$V4), split = ";")
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
setwd(refseq.path)
exons_region_refseq <- as_tibble(read.table("refseq_exons_171007.bed", sep = "\t")) %>%
  select(V1, V2, V3) %>%
  rename(chromosome = V1, start = V2, end = V3) %>%
  distinct()
exons_region_ext <- exons_region_refseq %>%
  select(chromosome, start, end) %>%
  mutate(start = start - half_exon, end = end + half_exon) %>%
  distinct()

tss_refseq <- as_tibble(read.table("refseq_TSS_hg19_170929.bed",
                                   header = F, sep = '\t')) %>%
  select(V1, V2, V3, V4) %>%
  rename(chromosome = V1, start = V2, end = V3, tss = V4) %>%
  arrange(chromosome, start, end, tss) %>%
  distinct()
tss_region_ext <- tss_refseq %>%
  select(chromosome, start, end) %>%
  mutate(start = start - half_tss, end = end + half_tss) %>%
  distinct()

data_enhancers_exclude <- data_enhancers_per %>% 
  genome_inner_join(tss_region_ext) %>%
  select(chromosome.x, start.x, end.x) %>%
  bind_rows(data_enhancers_per %>% 
              genome_inner_join(exons_region_ext) %>%
              select(chromosome.x, start.x, end.x)) %>%
  distinct() %>%
  rename(chromosome = chromosome.x, start = start.x, end = end.x)

data_enhancers_new <- data_enhancers_per %>%
  anti_join(data_enhancers_exclude)
e_all_new <- unique(data_enhancers_new$enhancer)
list_drop <- which(!e_eplinks_refseq %in% e_all_new) #691 rows to exclude

# Update the enhancer, tss, and gene lists
refseq_info_split_filtered <- refseq_info_split[-list_drop]
e_eplinks_refseq <- unlist(lapply(refseq_info_split_filtered, `[[`, 1)) 
tss_eplinks_refseq <- unlist(lapply(refseq_info_split_filtered, `[[`, 2))
g_eplinks_refseq <- unlist(lapply(refseq_info_split_filtered, `[[`, 3))
tss_all <- unique(unlist(strsplit(tss_eplinks_refseq, split = ",")))

# Build tss-enhancer links
pelinks_raw <- tibble(enhancer = e_eplinks_refseq,
                      tss = tss_eplinks_refseq, 
                      gene = g_eplinks_refseq
)

mapping_tss_enhancer <- pelinks_raw %>% 
  select(tss, enhancer) %>%
  group_by(tss) %>%
  summarise(enhancer = paste0(enhancer, collapse=";"))

mapping_tss_gene <- pelinks_raw %>% 
  select(tss, gene) %>%
  distinct()

pelinks_raw <- mapping_tss_enhancer %>%
  inner_join(mapping_tss_gene, by = "tss") %>%
  distinct() %>%
  select(-gene)
# these genes are annotated by fantom file enhancer_tss_associations.bed, drop for now

# Locate the TSSs
sum(!tss_all %in% tss_refseq$tss) #260 TSSs that no longer exist in the current RefSeq database
pelinks_raw_unnest <- pelinks_raw %>%
  mutate(tss = strsplit(tss, ",")) %>%
  unnest(tss)

pelinks_sameTSS <- tss_refseq %>% 
  left_join(pelinks_raw_unnest, by = "tss") %>%
  group_by(chromosome, start, end, enhancer) %>%
  summarise(tss = paste0(tss, collapse = ";")) %>%
  ungroup() %>%
  mutate(e_count = ifelse(is.na(enhancer), 0, sapply(strsplit(enhancer, ";"), length)))

# A histogram showing the number of enhancers each tss is associated with
setwd(output.path)
pelinks_plot <- pelinks_sameTSS %>%
  mutate(category = cut(e_count, breaks = c(0, 5, 10, 15, 20, Inf), right = F)) %>%
  group_by(category) %>%
  summarise(count = n())
png(file = "ep_hist.png", width = 1200, height = 800, res = 160)
ggplot(pelinks_plot, aes(x = category, y = count)) +
  geom_bar(stat = "identity", fill = "lightblue", alpha = 0.7, width = 0.5) + 
  geom_text(aes(label = count), nudge_y = 500) +
  labs(x = "number of enhancers per tss") +
  theme_bw()
dev.off()


## Part 3. Define promoters as the regions close to TSS but excluding exons
# check if exons overlap
exons_check_result <- nrow(exons_region_refseq %>% 
                             genome_inner_join(exons_region_refseq) %>% 
                             filter(start.x != start.y | end.x != end.y) %>% 
                             filter(start.x != end.y & start.y != end.x)) # due to bed file interval format ( , ]
print(ifelse(exons_check_result, 
             "there is overlap between exons", 
             "no overlap between exons")) # so there is overlap between exons, and we want to merge the overlapped ones
exons_region_refseq_me <- NULL # mutually exclusive exon regions
exons_update <- exons_region_refseq
while(nrow(exons_update) != 0){
  exons_self_join <- exons_update %>%
    genome_inner_join(exons_update)
  
  exons_no_overlap <- exons_self_join %>%
    group_by(chromosome.x, start.x, end.x) %>%
    count() %>%
    filter(n == 1) %>%
    select(-n) %>%
    rename(chromosome = chromosome.x,
           start = start.x,
           end = end.x)
  
  exons_overlap <- exons_self_join %>%
    filter(start.x < start.y)
  
  exons_update <- exons_overlap %>%
    mutate(end.x = pmax(end.x, end.y)) %>%
    select(1:3) %>%
    rename(chromosome = chromosome.x,
           start = start.x,
           end = end.x) %>%
    distinct()
  
  exons_region_refseq_me <- bind_rows(exons_region_refseq_me, exons_no_overlap)
}

data_promoters_raw <- tibble(chromosome = pelinks_sameTSS$chromosome, 
                             start = pelinks_sameTSS$start - half_p, 
                             end = pelinks_sameTSS$start + half_p,
                             tss = pelinks_sameTSS$tss)

promoters_join_exons <- data_promoters_raw %>%
  genome_left_join(exons_region_refseq_me)

promoters_without_exons <- promoters_join_exons %>%
  filter(is.na(chromosome.y)) %>%
  select(1:4) %>%
  rename(chromosome = chromosome.x,
         start = start.x,
         end = end.x)

promoters_with_exons <- promoters_join_exons %>%
  na.omit() %>%
  arrange(chromosome.x, start.x, end.x, start.y)

temp_regions <- promoters_with_exons %>%
  select(chromosome = chromosome.x, 
         start = start.x, 
         end = end.x,
         tss) %>%
         mutate(label = 1) %>%
  bind_rows(promoters_with_exons %>%
              select(chromosome = chromosome.y, 
                     start = start.y, 
                     end = end.y,
                     tss) %>%
                     mutate(label = 2)) %>%
  distinct() %>%
  arrange(chromosome, start, end)

data_promoters_noexons <- promoters_without_exons
for(tss_select in unique(temp_regions$tss)){
  temp <- temp_regions %>%
    filter(tss == tss_select)
  start_new <- c(temp$start[1], temp$end[-1])
  end_new <- c(temp$start[-1], temp$end[1])
  data_append <- tibble(
    chromosome = temp$chromosome,
    start = start_new,
    end = end_new,
    tss = temp$tss
  ) %>% filter(start < end)
  data_promoters_noexons <- bind_rows(data_promoters_noexons, data_append)
}

# Add the promoter variable in pelinks_sameTSS
pelinks_sameTSS_complete <- pelinks_sameTSS %>% 
  inner_join(data_promoters_noexons, by = "tss") %>%
  select(chromosome = chromosome.x, tss_pos = start.x, enhancer, tss,
         e_count, start = start.y, end = end.y) %>%
  mutate(promoter = paste0(chromosome, ":", start, "-", end)) %>%
  group_by(chromosome, tss, enhancer, e_count, tss_pos) %>%
  summarise(promoter = paste0(promoter, collapse = ";")) %>%
  ungroup() %>%
  select(tss, chromosome, tss_pos, promoter, enhancer, e_count)

promoters_associated <- pelinks_sameTSS_complete %>% select(tss, promoter)
enhancers_associated <- pelinks_sameTSS_complete %>% select(tss, enhancer)


# ## Part 4. Annotate the coding exons for each TSS
# setwd(refseq.path)
# refseq_codingexons <- 
#   as_tibble(read.table("refseq_codingexons_171007.bed", sep = "\t")) %>%
#   select(V1, V2, V3, V4) %>%
#   rename(chromosome = V1, start = V2, end = V3, tss = V4) %>%
#   arrange(chromosome, start, end, tss) %>%
#   distinct() %>%
#   mutate(tss = gsub("_cds.*$", "", tss))

## Output the result: TSS, promoter, enhancer, coding exons, etc.
setwd(output.path)
write.table(pelinks_sameTSS_complete, "pelinks.txt", quote = F, row.names = F)
write.csv(pelinks_sameTSS_complete, "pelinks.csv", quote = F, row.names = F)

