library(tidyverse)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(fuzzyjoin)

## arguments from the bash file ##
args <- commandArgs(trailingOnly = TRUE)
cmd_arg <- as.numeric(args[1])

# ============ PATH INFORMATION ============ 
setwd("/Users/vincentiusmartin/Research/CancerMutation")
source("functions.R")

promoter_path <- "input/promoters_refseq_coding.txt"
enhancer_path <- "input/all_enhancers_fantom.txt"
exon_path <- "input/all_exons_refseq.txt"
icgc_snps_path <- "input/icgc_snps_LIRI_JP.tsv"

outdir <- "output"
cistype <- "enhancer" # "enhancer" or "promoter"

# ============ READ INPUT FILES ============ 
icgc_snps <- fread(icgc_snps_path)

promoters <- fread(promoter_path, col.names = c("chromosome", "start", "end", "promoter","strand")) %>%
  select(chromosome, start, end, promoter) %>% distinct()
enhancers <- fread(enhancer_path, col.names = c("chromosome", "start", "end", "enhancer"))
trimut_df <- fread(paste0("input/",cistype,"_mut_rate.csv"), select=c("ref_tri","mut_tri","tri_mut_rate"))
# ============ PROCESSING ============ 
cis <- if (cistype == "enhancer") enhancers else promoters
genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg19")

# get genomic regions with mutations
cis_snps_join <- icgc_snps %>% 
  genome_inner_join(cis, by = c("chromosome", "start", "end")) 

# this contains cis with snps
cis_snps <- cis_snps_join %>%
  select(chromosome = chromosome.y,
         start = start.y,
         end = end.y,
         !!cistype) %>%
  arrange(chromosome,start,end) %>%
  distinct()

# all snps that happened within the cis-regions
snps_in_cis <- cis_snps_join %>%
  select(chromosome = chromosome.x,
         pos = start.x,
         cis_chromosome = chromosome.y,
         cis_start = start.y,
         cis_end = end.y,
         ref, mut) %>%
  distinct()

# the cis bg sequence including 5bp context
cis_seqs_bg <- getSeq(genome, names = cis_snps$chromosome,
                      start = cis_snps$start - 5,
                      end = cis_snps$end + 5)

result_df <- lapply(seq(nrow(cis_snps)), function(i){
  ## Step 1: Generate all possible single mutations for enhancers with mutations ##
  # get the location and sequence of the enhancer
  cur_cis <- cis_snps[i,]
  cur_cis_seq <- subseq(cis_seqs_bg[[i]], start = 6, end = -6) # we take the actual cis-sequence from the bg
  cur_cis_11mer <- DNAStringSet(cis_seqs_bg[[i]], 
                                start = 1:(length(cis_seqs_bg[[i]])-10),  
                                width = 11) # cut the cis into 11mers
  
  # generate all possible single mutations
  cis_all_possible_mutations <- tibble(
    chromosome = cur_cis$chromosome,
    pos = seq(cur_cis$start, cur_cis$end),
    seq_mut_bg = as.character(cur_cis_11mer)
  ) %>% 
    mutate(ref = unlist(strsplit(as.character(cur_cis_seq), ""))) %>% 
    dplyr::slice(rep(seq(length(cur_cis_seq)), each = 4)) %>% 
    mutate(mut = rep(c("A", "C", "G", "T"), length(cur_cis_seq))) %>% 
    filter(ref != mut) %>% 
    mutate(twimer = paste0(seq_mut_bg, mut)) %>%
    rowwise() %>%
    mutate(idx = seqtoi(twimer)) %>%
    ungroup() %>% 
    select(-seq_mut_bg)
  
  # mut_rate is  normalized so that they add up to 1.
  cis_snps_bg_prediction <- cis_all_possible_mutations %>%
    mutate(ref_tri = substr(twimer, 5, 7),
           mut_tri = paste0(substr(ref_tri, 1, 1), mut, substr(ref_tri, 3, 3)),
           ref_tri_rev = as.character(reverseComplement(DNAStringSet(ref_tri))),
           mut_tri_rev = as.character(reverseComplement(DNAStringSet(mut_tri))),
           ref_tri = ifelse(ref_tri < ref_tri_rev, ref_tri, ref_tri_rev),
           mut_tri = ifelse(ref_tri < ref_tri_rev, mut_tri, mut_tri_rev)) %>% 
    select(-c(ref_tri_rev, mut_tri_rev))  %>% 
    inner_join(trimut_df, by = c("ref_tri", "mut_tri")) %>% 
    mutate(cond_tri_mut_rate = tri_mut_rate/sum(tri_mut_rate))
  
  # mutation effect prediction for the actual mutations in the cis
  mut_effect_actual_snps <- snps_in_cis %>% 
    filter(cis_chromosome==cur_cis$chromosome, 
           cis_start==cur_cis$start, 
           cis_end==cur_cis$end) %>%
    select(chromosome, pos, ref, mut) %>%
    inner_join(cis_snps_bg_prediction, by = c("chromosome", "pos", "ref", "mut"))
  
  ## Step 3: get the cis level mutation effect and p-value ##
  cis_diff <- absmax(mut_effect_actual_snps$diff)
  n_snps <- nrow(mut_effect_actual_snps)
  
  # a subset of all possible cis mutations 
  # with abs(mutation effect) < abs(mut_effect_enhancer)
  diff_less_than_actual <- cis_snps_bg_prediction %>%
    filter(abs(diff) < abs(cis_diff))
  
  # two-sided p-value of seeing more extreme cis level mutations effect
  ts_p_value <- 1 - sum(diff_less_than_actual$cond_tri_mut_rate)^n_snps
  
  # approximate one-sided p-value of seeing cis effect < mut_effect_cis
  os_p_value <- ifelse(cis_diff < 0, ts_p_value/2, 1 - ts_p_value/2)
  
  result <- c(cis = cur_cis[[cistype]],
              mut_effect = cis_diff,
              p_less = os_p_value) %>% rename(cis = cistype)
  return(result)                    
}) %>% bind_rows