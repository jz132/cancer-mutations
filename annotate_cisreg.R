library(tidyverse)
library(data.table)

setwd("/Users/vincentiusmartin/Research/CancerMutation")

cistype <- "promoter" 
snp_count_dir <- paste0("dataset/epdata/",cistype,"_snp_count.csv")

# for enhancers
eplink_path <- "dataset/input/eplinks-fantom-filtered.csv" # enhancer-promoter to gene link

# for promoters
hgnc_map_path <- "dataset/input/RefSeq/refseq_to_hgnc.txt"
refseq_tss_map_path <- "dataset/input/RefSeq/refseq_TSS_hg19_170929.bed"

# ============ Map the SNPs data ============ 
snp_count <- fread(snp_count_dir) %>% filter(count > 0) %>% select(-length)

if(cistype == "promoter"){
  refseq_tss_mapping <- fread(refseq_tss_map_path, select=c("V1","V2","V4"), col.names=c("chromosome","pos","tss")) %>%
    filter(grepl("NM", tss))
  
  hgnc_mapping <- fread(hgnc_map_path) %>%
    select(gene = "Approved symbol",
           refseq_id = "RefSeq IDs",
           refseq_id_ncbi = "RefSeq(supplied by NCBI)") %>%
    pivot_longer(cols = refseq_id:refseq_id_ncbi, values_to = "tss") %>%
    select(gene, tss) %>%
    filter(grepl("NM", tss)) %>%
    na.omit() %>%
    distinct()
  
  # ============ Annotate SNPs with the gene names ============ 
  snp_annotated <- snp_count %>%
    mutate(pos = end - 1000) %>%
    inner_join(refseq_tss_mapping, by = c("chromosome", "pos")) %>%
    inner_join(hgnc_mapping, by = "tss") %>%
    select(-c("pos", "tss")) %>%
    distinct() 
}else{ # enhancer
  eplinks <- fread_sep(eplink_path, "enhancer") %>% distinct(chromosome,start,end,gene)
  snp_annotated <- snp_count %>% 
    inner_join(eplinks, by=c("chromosome", "start", "end"))
}

write.csv(snp_annotated, paste0(cistype,"_annotated.csv"), quote=FALSE, row.names=FALSE)