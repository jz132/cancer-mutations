library(data.table)
library(tidyverse)

# ============ PATH INFORMATION ============ 
source("functions.R")

promoter_path <- "dataset/epdata/all_promoters_refseq.txt"
refseq_tss_path <- "dataset/epdata/refseq_TSS_hg19_170929.bed"
hgnc_map_path <- "dataset/epdata/refseq_to_hgnc.txt"
outpath <- "dataset/epdata/promoters_refseq_coding.txt"

# ============ PROCESS ============ 
promoters <- fread(promoter_path, col.names = c("chromosome", "start", "end", "promoter", "strand")) %>%
  distinct()

# keep only promoters for protein coding genes
mapping_refseq_tss <- fread(refseq_tss_path) %>%
  select(c(1,2,4)) %>%
  dplyr::rename(chromosome = "V1", pos = "V2", tss = "V4") %>% 
  filter(grepl("NM", tss))

mapping_hgnc <- fread(hgnc_map_path) %>%
  select(gene = "Approved symbol",
         refseq_id = "RefSeq IDs",
         refseq_id_ncbi = "RefSeq(supplied by NCBI)") %>%
  pivot_longer(cols = refseq_id:refseq_id_ncbi, values_to = "tss") %>%
  select(gene, tss) %>%
  filter(grepl("NM", tss)) %>%
  na.omit() %>%
  distinct()

promoters_coding <- promoters %>% 
  mutate(pos = (start + end)/2) %>% 
  inner_join(mapping_refseq_tss) %>%
  inner_join(mapping_hgnc) %>%
  select(-c("pos", "tss", "gene")) %>% 
  distinct()

nrow(promoters_coding %>% select(chromosome,start,end) %>% distinct())
write.table(promoters_coding, outpath, quote = F, row.names = F, col.names = F, sep = "\t")
