library(tidyverse)
library(data.table)
options(dplyr.summarise.inform = FALSE)
rename <- dplyr::rename

setwd("/Users/vincentiusmartin/Research/CancerMutation")
source("functions.R")

# for 1 regulatory region, just use vector of length 1
tf_mapping_path <- "dataset/input/QBiC_mapping_582_TFgenes_to_666_uPBM_datasets.csv"
pred_dirs <- c("dataset/predictions/preds_lirijp_enhancers", "dataset/predictions/preds_lirijp_promoters")  #c("dataset/predictions/preds_lirijp_enhancers", "dataset/predictions/preds_lirijp_promoters")
entropy_paths <- c("dataset/epdata/enhancer_entropy.csv", "dataset/epdata/promoter_entropy.csv") #c("dataset/epdata/enhancer_entropy.csv", "dataset/epdata/promoter_entropy.csv")
cis_annotated_paths <- c("dataset/annotated/enhancer_annotated.csv", "dataset/annotated/promoter_annotated.csv") #c("dataset/annotated/enhancer_annotated.csv", "dataset/annotated/promoter_annotated.csv")
cistypes <- c("enhancer", "promoter") #c("enhancer", "promoter")
cislen <- length(cistypes)

# ST3GAL6
# ============ Read SNPs data ============ 
cis_annotated <- lapply(cis_annotated_paths,fread) 

cis_annotated_wtype <- cis_annotated
for(i in seq(cislen)){
  cis_annotated_wtype[[i]]$type <- cistypes[i]
}
cis_annotated_combined <- bind_rows(cis_annotated_wtype)

# ============ Annotate SNPs with the gene names ============ 
genes_count <- cis_annotated_combined %>%
  select(gene,count) %>%
  group_by(gene) %>%
  summarise(count = sum(count)) %>%
  arrange(gene)
gene_names <- genes_count %>% select(gene)

# ============ Entropy data ============ 
entropy_df <- mapply(fread_sep, entropy_paths, cistypes, SIMPLIFY = FALSE) 
cis_df <- lapply(seq(cislen), function(i){
  return (cis_annotated[[i]] %>% 
    inner_join(entropy_df[[i]], by=c("chromosome","start","end"))
)})

# ============ Mapping data ============ 
tf_mapping <- fread(tf_mapping_path, select=c("gene","upbm")) %>%
  mutate(upbm = strsplit(upbm,";")) %>%
  unnest(cols=c(upbm)) %>%
  mutate(upbm = tools::file_path_sans_ext(trimws(upbm))) %>%
  rename(tf = gene) %>%
  group_by(upbm) %>%
  summarise(tf = paste0(tf, collapse = ";"))

# ============ Directory with predictions ============ 
# we can just use one of the prediction directory as they have the same pbmnames
allpreds <- lapply(pred_dirs, function(path){Sys.glob(paste0(path,"/*.csv"))})
end_pattern_pos <- as.integer(str_locate(allpreds[[1]][1],"prediction6mer.")[,2]) + 1
pbmnames <- tools::file_path_sans_ext(str_sub(allpreds[[1]], start=end_pattern_pos))

#which(str_detect(allpreds[[1]], "RUNX1"))
# iterate over the pbm files
maxiter <- length(allpreds[[1]])
result_zscores <- lapply(seq(maxiter), function(i){
  
  # might need to put mean and sd to qnorm
  mutated_effects_cis <- lapply(seq(cislen), function(cis_idx){
    return(fread_sep(allpreds[[cis_idx]][i], cistypes[cis_idx]) %>%
      inner_join(cis_df[[cis_idx]], by = c("chromosome", "start", "end")) %>%
      mutate(z_score = qnorm(p_less), cistype = cistypes[cis_idx]) %>% 
      select(-p_less)
  )}) %>% bind_rows()
  
  ## for example analyses only
  # cis_df[[1]] %>% filter(start == 230177949)
  # mutated_effect_summarised <- mutated_effects_cis %>% 
  #   group_by(gene) %>%
  #   summarise(count = sum(count),
  #             row_count = n(),
  #             z_score = sum(entropy*z_score)/sqrt(sum(entropy^2))) %>%
  #   arrange(desc(z_score))
  # mutated_effect_summarised %>% filter(gene == "ATXN3")
  # mutated_effects_cis %>% filter(gene == "ATXN3")
  ##
  
  # Combine p-values using Stouffer's Z-score method
  mutated_effect_summarised <- mutated_effects_cis %>% 
    group_by(gene) %>%
    summarise(count = sum(count),
              z_score = sum(entropy*z_score)/sqrt(sum(entropy^2))) %>%
    arrange(gene)
  # to make sure we have the correct ordering
  # mutated_effect_summarised <- gene_names %>% inner_join(mutated_effect_summarised, by="gene")
  if(nrow(mutated_effect_summarised) != nrow(gene_names)){
    print(paste("Different gene count for",allpreds[i]))
  }
    
  return (list(z_cis = mutated_effects_cis$z_score, z_gene = mutated_effect_summarised$z_score))
})

# ======== CIS LEVEL ======== 
z_cis_all <- lapply(result_zscores[1:(length(result_zscores))], '[[' , "z_cis") 
names(z_cis_all) <- pbmnames[1:(length(pbmnames))]
z_cis_df <- bind_cols(z_cis_all)
p_cis_df <- z_cis_df %>% mutate_all(~ 2*(1-pnorm(abs(.)))) # 2 for two sided p-value
p_cis_adjusted_df <- apply(p_cis_df, 2, p.adjust, method = "hochberg")

# cis level results
z_score_absmax_per_cis <- apply(z_cis_df, 1, absmax)
p_min_per_cis <- apply(p_cis_df, 1, min) #each row is a cis-region
adjusted_p_min_per_enhancer <- apply(p_cis_adjusted_df, 1, min) #each row is a cis-region
upbm_min_per_cis <- pbmnames[apply(p_cis_df, 1, which.min)]
mutated_max_effect_per_cis <- cis_annotated_combined %>% 
  mutate(absmax_z_score = z_score_absmax_per_cis,
         min_p_value = p_min_per_cis,
         min_adjusted_p_value = adjusted_p_min_per_enhancer,
         upbm = upbm_min_per_cis) %>%
  inner_join(tf_mapping, by="upbm") %>%
  arrange(min_p_value) 
hist(log10(mutated_max_effect_per_cis$min_p_value))
summary(log10(mutated_max_effect_per_cis$min_p_value))

cisname <- ifelse(length(cistypes) == 1, cistypes[1], "combined")
write.csv(mutated_max_effect_per_cis, paste0(cisname,"_muteffect.csv"), quote=FALSE, row.names=FALSE)

# ======== GENE LEVEL ======== 
z_genes_all <- lapply(result_zscores, '[[' , "z_gene") 
names(z_genes_all) <- pbmnames
z_genes_df <- bind_cols(z_genes_all)
p_genes_df <- z_genes_df %>% mutate_all(~ 2*(1-pnorm(abs(.))))
p_genes_adjusted_df <- apply(p_genes_df, 2, p.adjust, method = "hochberg")

hist(p_genes_df$`Homo_sapiens|NA|Shen2018|MYC-MAX`)

# save the z-score tables
z_mat_all <- z_genes_df %>% 
  mutate(gene = genes_count$gene) %>% 
  select(gene, everything()) 
# write.csv(format(z_mat_all,digits = 2), paste0(cisname,"_zscore_allgenes.csv"), quote=FALSE, row.names=FALSE)

# gene level results 
p_min_per_gene <- apply(p_genes_df, 1, min)
p_min_per_gene_adjusted <- apply(p_genes_adjusted_df, 1, min)
upbm_min_per_gene <- pbmnames[apply(p_genes_df, 1, which.min)]

max_effect_per_gene <- genes_count %>%
  mutate(min_p_value = p_min_per_gene,
         min_adjusted_p_value = p_min_per_gene_adjusted,
         upbm = upbm_min_per_gene) %>% 
  inner_join(tf_mapping, by="upbm") %>%
  arrange(min_p_value)
hist(log10(max_effect_per_gene$min_p_value))
summary(log10(max_effect_per_gene$min_p_value))

write.csv(max_effect_per_gene, paste0(cisname,"_gene_muteffect.csv"), quote=FALSE, row.names=FALSE)
