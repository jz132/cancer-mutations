rename <- dplyr::rename
library(fuzzyjoin)
source("functions.R")

library(TFBSTools)
suppressMessages(library(JASPAR2020))

library(BSgenome.Hsapiens.UCSC.hg19)
genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg19")

reg_muteff_path <-  "dataset/combined_predictions/combined_muteffect.csv"
icgc_snps_path <- "dataset/input/icgc_snps_LIRI_JP_no_exons.tsv"
pred_dir <- "dataset/predictions"

genes_target <- c("C1R","CENPA","DPM3","ATXN3","CEBPD","ZNF561","TM4SF18","CTNNA3")
# genes_target <- c("RNF32","DDX21","N4BP2L2","ZNF606","ST3GAL6","ACBD5","NLRP3","JAK1",
#           "TET2","RNF39","CXCR2","CNTD1","UBAP2L",
#           "SYT1","TBC1D14","PPM1M","C3orf14","ARL4C","TSNAXIP1","B4GALNT2","TAF6","SF1","SPACA3","NLRC4")

# ======== Read input files ========   

reg_muteff <- fread(reg_muteff_path) %>% 
  filter(gene %in% genes_target) %>%
  select(chromosome, start, end, gene, upbm, tf, type) %>%
  rename(cis_type = type)

icgc_snps <- fread(icgc_snps_path) %>% select(chromosome, start, end, ref, mut, icgc_donor_id)

# ======== Processing ========   

top_cis_with_snps <- reg_muteff %>% 
  genome_inner_join(icgc_snps, by=c("chromosome", "start", "end")) %>%
  rename( chromosome = "chromosome.x",
          start = "start.x",
          end = "end.x",
          snp_pos = "start.y"
  ) %>%
  select(-c("chromosome.y","end.y")) %>%
  arrange(gene) %>%
  distinct()

top_cis_diff <- lapply(seq(nrow(top_cis_with_snps)), function(i){
  cur_cis <- top_cis_with_snps[i,]
  cur_pred_path <- Sys.glob(paste0(pred_dir,"/preds_lirijp_",cur_cis$cis_type,"s/*",cur_cis$upbm,"*"))
  cur_pred <- fread_sep(cur_pred_path,cur_cis$cis_type)
  return(cur_cis %>% inner_join(cur_pred, by=c("chromosome","start","end")))
})  %>% bind_rows()

windowlen <- 5
seqs_top_cis <- as.character(
  getSeq(genome, names = top_cis_diff$chromosome,
         start = top_cis_diff$snp_pos - windowlen,
         end = top_cis_diff$snp_pos + windowlen)
)

top_cis_diff$seq <- seqs_top_cis
write.csv(top_cis_diff %>% distinct(), paste0("top_cis_diff.csv"), quote=FALSE, row.names=FALSE)

# how many promoter only genes have enhancer mutations
# length(unique((top_cis_diff %>% select(gene,cis_type) %>% filter(cis_type=="enhancer"))$gene))

tcdiff_split <- top_cis_diff %>%
  mutate(tf = strsplit(tf,";")) %>%
  unnest(cols=c(tf)) %>%
  distinct()

# allsites <- lapply(seq(nrow(tcdiff_split)), function(i){
#   curtf <- tcdiff_split[i,]$tf
#   writeLines(paste(i,curtf))
#   curseq <- tcdiff_split[i,]$seq
#   snp_pos <- tcdiff_split[i,]$snp_pos - tcdiff_split[i,]$start
#   
#   opts <- list("species"=9606, "name"=curtf)
#   pfm <- getMatrixSet(JASPAR2020, opts)
#   if(length(pfm) == 0){
#     pfm <- getMatrixSet(JASPAR2020, list("name"=curtf))
#     if(length(pfm) == 0){
#       return(data.frame(list(tf=curtf, seq=curseq, site_start=-999,site_end=-999,absScore=-999)))
#     }
#   }
#   pfm <- pfm[[1]]
#   
#   res <- as.data.frame(searchSeq(toPWM(pfm), curseq, min.score = "80%"))
#   
#   if(nrow(res) == 0){
#     return(data.frame(list(tf=curtf, seq=curseq, site_start=-999,site_end=-999,absScore=0)))
#   }else{
#     res$tf = curtf
#     res$seq=curseq
#     res <- res %>% rename(site_start = start, site_end = end) %>%
#       select(tf,seq,site_start, site_end, absScore) %>% distinct()
#     return(res)
#   }
# }) %>% bind_rows()
# 
# tcdiff_split <- tcdiff_split %>% 
#   inner_join(allsites, by=(c("tf","seq"))) %>% 
#   rename(pwm_score="absScore") %>%
#   distinct()
# write.csv(tcdiff_split, paste0("top_cis_diff_wsites.csv"), quote=FALSE, row.names=FALSE)

# ======= Check the binding on specific TFs # ======= 
gene_target <- "ATXN3"

filename <- "Homo_sapiens|NA|Shen2018|RUNX1.csv"
tf_enhancer_preds <- fread_sep(paste0(
               "dataset/predictions/preds_lirijp_enhancers/enhancers_mutation_result_LIRI-JP_prediction6mer.", 
               filename),"enhancer")
tf_enhancer_preds$cis_type <- "enhancer"
tf_promoter_preds <- fread_sep(paste0(
  "dataset/predictions/preds_lirijp_promoters/promoter_mutation_result_LIRI-JP_prediction6mer.", 
  filename),"promoter")
tf_promoter_preds$cis_type <- "promoter"

all_reg_preds <- bind_rows(tf_enhancer_preds, tf_promoter_preds)

geneseqs <- tcdiff_split %>% 
  filter(gene == gene_target) %>%
  select(chromosome,start,end,cis_type,icgc_donor_id,seq,mut) %>%
  distinct() %>% inner_join(all_reg_preds) %>%
  select(chromosome,start,end,cis_type,icgc_donor_id,mut_effect,p_less)
