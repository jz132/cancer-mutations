library(dplyr)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(fuzzyjoin)

slice <- dplyr::slice

## arguments from the bash file ##
args <- commandArgs(trailingOnly = TRUE)
cmd_arg <- 1 #as.numeric(args[1])

setwd("/Users/vincentiusmartin/Research/CancerMutation")
source("functions.R")
# ============ PATH INFORMATION ============ 
icgc_snps_path <- "dataset/input/icgc_snps_LIRI_JP_no_exons.tsv"
promoter_path <- "dataset/input/promoters_refseq_coding.txt"
enhancer_path <- "dataset/input/all_enhancers_fantom.txt"
cistype <- "enhancer" # "enhancer" or "promoter"

tbl12mer_dir <- "pbm-all"
tbl12mer_path <- Sys.glob(paste0(tbl12mer_dir,"/*.txt"))[cmd_arg]
cat(paste0("Processing ",tbl12mer_path,"\n"))
tbl12mer <- fread(tbl12mer_path, fill=TRUE) 

# ============ READ INPUT FILES ============ 
icgc_snps <- fread(icgc_snps_path)
num_donors <- nrow(icgc_snps %>% distinct(icgc_donor_id))

cis <- if (cistype == "enhancer"){
  fread(enhancer_path, col.names = c("chromosome", "start", "end", "enhancer"))
}else{
  fread(paste0("dataset/epdata/",cistype,"_mut_rate.csv"), select=c("ref_tri","mut_tri","tri_mut_rate"))
}

trimut_df <- fread(paste0("dataset/epdata/",cistype,"_mut_rate.csv")) 
ref_mut_idx <- trimut_df %>% select(ref_tri,mut_tri) %>% mutate(ref_mut_idx = row_number())

mut_rate_df <- trimut_df %>% select(ref_tri,mut_rate) %>% distinct()
trimutrate <- mut_rate_df$mut_rate
names(trimutrate) <- mut_rate_df$ref_tri

# ============ PROCESSING ============ 

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
         cis_start = start.y,
         cis_end = end.y,
         ref, mut) %>%
  distinct()

###

genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg19")
cis_seqs_bg <- getSeq(genome, names = cis_snps$chromosome,
                      start = cis_snps$start - 5,
                      end = cis_snps$end + 5)

# simulation Setting
set.seed(1234)
n_run <- 10^2

div_const <- nrow(cis_snps) %/% 100 # for progress report
result_df <- lapply(seq(nrow(cis_snps)), function(i){
  if (i %% div_const == 0) {
    writeLines(paste0("Progress ",i,"/",nrow(cis_snps)))
  }
    
  cur_cis <- cis_snps[i,]
  cur_cis_seq <- subseq(cis_seqs_bg[[i]], start = 6, end = -6) # we take the actual cis-sequence from the bg
  cur_cis_11mer <- DNAStringSet(cis_seqs_bg[[i]], 
                                start = 1:(length(cis_seqs_bg[[i]])-10),  
                                width = 11) # cut the cis into 11mers
  
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
  
  cis_snps_bg_prediction <- cis_all_possible_mutations %>% 
    mutate(ref_tri = substr(twimer, 5, 7),
           mut_tri = paste0(substr(ref_tri, 1, 1), mut, substr(ref_tri, 3, 3)),
           ref_tri_rev = as.character(reverseComplement(DNAStringSet(ref_tri))),
           mut_tri_rev = as.character(reverseComplement(DNAStringSet(mut_tri))),
           ref_tri = ifelse(ref_tri < ref_tri_rev, ref_tri, ref_tri_rev),
           mut_tri = ifelse(ref_tri < ref_tri_rev, mut_tri, mut_tri_rev),
           diff = tbl12mer$diff[idx]
           ) %>% 
    select(-c(ref_tri_rev, mut_tri_rev)) %>%
    inner_join(ref_mut_idx, by = c("ref_tri", "mut_tri"))
  
  ## Simulation
  
  # the count of each trinucleotide in the enhancer
  vec_tri <- reverseMerge_vec(trinucleotideFrequency(cur_cis_seq))
  n_tri <- length(vec_tri)
  
  #  === first random sampling from a binomial distribution for the number of mutated trinucleotide ===
  # Find 'n_run' random values from a sample of 'vec_tri[j]*num_donors' with probability of 'trimut_rate
  # Get how many trinucleotides are mutated on each run
  ref_sim <- lapply(names(vec_tri), function(curtri){
    rbinom(n_run, vec_tri[curtri]*num_donors, trimutrate[curtri])  # no need to multiply if just for 1 donor
  }) 
  names(ref_sim) <- mut_rate_df$ref_tri
  ref_sim <- ref_sim %>% bind_cols()
  
  rbinom(4,100,0.3)
  
  # === then random sampling from a multinomial distribution for each trinucleotide to trinucleotide mutation ===
  
  # probability of mutating to another base
  mut_prop_df <- matrix(trimut_df %>%
    select(ref_tri, mut_tri, proportion) %>% 
    distinct() %>%
    pull(proportion), nrow = 32, ncol = 3, byrow = T)
  rownames(mut_prop_df) <- mut_rate_df$ref_tri
  
  ref_mut_sim <- lapply(mut_rate_df$ref_tri, function(tri){
    return(t(sapply(ref_sim[[tri]], function(n){
      return(t(rmultinom(1, n, mut_prop_df[tri,])))
    })))
  })
  ref_mut_sim <- as.data.frame(ref_mut_sim)
  colnames(ref_mut_sim) <- 1:ncol(ref_mut_sim)
  
  ref_mut_sim_diff <- lapply(ref_mut_idx$ref_mut_idx, function(rmidx){
    cur_diffs <- (cis_snps_bg_prediction %>% filter(ref_mut_idx == rmidx))$diff
    return(sapply(ref_mut_sim[,rmidx], function(x){
      return(ifelse(x != 0, absmax(sample(cur_diffs, x, replace = T)),0))
    }))
  })
  ref_mut_sim_diff <- do.call(cbind,ref_mut_sim_diff)
  
  diff_actual_snps <- snps_in_cis %>% 
    filter(chromosome==cur_cis$chromosome, 
           cis_start==cur_cis$start, 
           cis_end==cur_cis$end) %>%
    select(chromosome, pos, ref, mut) %>%
    inner_join(cis_snps_bg_prediction, by = c("chromosome", "pos", "ref", "mut"))
  
  cis_diff <- absmax(diff_actual_snps$diff)
  mut_enhancer_sim_effect <- apply(ref_mut_sim_diff, 1, absmax)
  p_value <- sum(abs(mut_enhancer_sim_effect) >= abs(cis_diff))/n_run # average probability across each run
  
  #print(paste("empirical p-value: ", p_value))
  result <- tibble({{cistype}} := cur_cis[[cistype]],
                   mut_effect = cis_diff,
                   p_less = p_value)
  return(result)
}) %>% bind_rows()
