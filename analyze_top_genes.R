library(tidyverse)
library(data.table)
library(fuzzyjoin)
source("functions.R")

rename <- dplyr::rename

setwd("/Users/vincentiusmartin/Research/CancerMutation")

# ======== Input path ========  

gene_muteff_all_path <- c("dataset/combined_predictions/promoter_gene_muteffect_expanded.csv",
                          "dataset/combined_predictions/enhancer_gene_muteffect_expanded.csv",
                          "dataset/combined_predictions/combined_gene_muteffect_expanded.csv")
regelms <- c("promoter","enhancer","combined")

pathology_path <- "dataset/validation/pathology.tsv"

# ======== Read input files ========  

# `prognostic - favorable`
pathology <- fread(pathology_path) %>% 
  filter(!is.na(`prognostic - favorable`) | !is.na(`prognostic - unfavorable`)) %>%
  select(`Gene name`,Cancer) %>% 
  `colnames<-`(c("gene", "cancer")) %>%
  group_by(gene) %>%
  mutate(cancer = paste0(cancer, collapse = ",")) %>%
  distinct()

gene_muteff <- lapply(seq(regelms), function(i){
  gm <- fread(gene_muteff_all_path[i]) %>% 
    arrange(min_adjusted_p_value) %>%
    mutate(neglog = -log(min_adjusted_p_value),
           rank = row_number())
  
  patho_genes <- gm %>% 
    filter(min_adjusted_p_value < 0.05) %>%
    select(gene) %>%
    left_join(pathology, by = "gene")
  
  return(list("gene_muteff" = gm, "patho_genes" = patho_genes))
})
names(gene_muteff) <- regelms

# ======== Fisher test of the pathogenic top genes ========   
for(curcis in regelms){
  cat(paste("Elm:",curcis,"\n"))
  
  cur_patho <- gene_muteff[[curcis]]$patho_genes 
  top_gene_patho_count <- nrow(cur_patho %>% filter(!is.na(cancer)))
  top_gene_patho_vec <- c(top_gene_patho_count,  nrow(cur_patho) - top_gene_patho_count)
  
  # use the median
  cur_gm <- gene_muteff[[curcis]]$gene_muteff
  gm_median <- median(cur_gm$min_p_value)
  gm_ctrl <- cur_gm %>% filter(min_p_value > gm_median) 
  ctrl_gene_patho_count <- length(intersect(gm_ctrl$gene, pathology$gene))
  ctrl_gene_patho_vec <- c(ctrl_gene_patho_count, length(gm_ctrl$gene) - ctrl_gene_patho_count)
  
  mat <- data.frame(
    "Top genes" = top_gene_patho_vec,
    "Ctrl genes" = ctrl_gene_patho_vec,
    row.names = c("Prognostic", "Non-prognostic"),
    stringsAsFactors = FALSE
  )
  
  print(mat)
  
  fisher_res <- fisher.test(mat, alternative = "greater")
  chisq_res <- chisq.test(mat, correct=FALSE)
  
  cat(paste("Fisher, p =",fisher_res$p,"\n"))
  cat(paste("Chisq, p =",chisq_res$p.value,"\n"))
  
  # Pathogenic for genes with significant mutated vs unmutated
  gm_top <- cur_gm %>% filter(min_adjusted_p_value < 0.05)
  gene_mut_unmut <- (gm_top %>% filter(p_value_mutated_vs_unmutated < 0.1))$gene
  gene_mut_patho <- intersect(gene_mut_unmut, pathology$gene)
  cat(paste("All top mut_unmut genes",paste(gene_mut_unmut, collapse = ','),"\n"))
  cat(paste("Significant pathogenic genes",paste(gene_mut_patho, collapse = ','),"\n"))
  
  gm_top <- gm_top %>%
    inner_join(cur_patho, by = c("gene")) %>%
    mutate(prognostic = ifelse(is.na(cancer),"prognostic","not prognostic")
    )
  gm_top$gene  <- factor(gm_top$gene, gm_top$gene)
  ggplot(gm_top, aes(x=gene, y=neglog, fill=prognostic)) + 
    geom_bar(stat="identity")+
    ylab("-log(p)") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    scale_fill_manual(values=c("#989898","#FFA500")) +
    theme_classic() +
    theme(axis.text.x=element_text(angle=45,hjust=1))
  ggsave(paste0(curcis,"_top_genes_barplot.pdf"), width = 10, height = 4)
  
  # p_value_mutated_vs_unmutated or p_value_cancer_vs_normal
  delta_exp_top <- -log(na.omit(gm_top$p_value_mutated_vs_unmutated))
  delta_exp_ctrl <- -log(na.omit(gm_ctrl$p_value_mutated_vs_unmutated))
  p_delta_exp <- wilcox.test(delta_exp_top,delta_exp_ctrl,"greater")
  cat(paste("P value patient mutated vs unmutated",p_delta_exp$p.value))
  
  b <- c(0, 0.05, 0.1, 1)  
  names <- c("p<=0.05", "0.05<p<=0.01", "p>0.01")
  cur_gm$pcat <- cut(cur_gm$min_adjusted_p_value, breaks = b, labels = names)
  
  gene_unmut_mut <- cur_gm %>% 
    select(p_value_mutated_vs_unmutated, pcat) %>% 
    filter(!is.na(p_value_mutated_vs_unmutated)) %>%
    mutate(neglog_unmut_mut = -log(p_value_mutated_vs_unmutated))
  
  ggplot(gene_unmut_mut, aes(x=pcat, y=neglog_unmut_mut, fill=pcat)) + 
    ylab("-log(p)") +
    geom_boxplot() +
    ylim(0,7.5) +
    scale_fill_manual(values=c("firebrick","mistyrose","#D3D3D3")) + 
    theme(legend.position = "none") 
  ggsave(paste0(curcis,"_exp_diff_unmut_mut.pdf"))
  
  # when hypothesis testing is needed
  g1 <- gene_unmut_mut %>% filter(pcat == "p<=0.05")
  g2 <- gene_unmut_mut %>% filter(pcat == "0.05<p<=0.01")
  g3 <- gene_unmut_mut %>% filter(pcat == "p>0.01")
  g1_v_g2 <- wilcox.test(g1$neglog_unmut_mut,g2$neglog_unmut_mut,"greater")
  g1_v_g3 <- wilcox.test(g1$neglog_unmut_mut,g3$neglog_unmut_mut,"greater")
  
  cat(paste("p<=0.05 vs 0.05<p<=0.01, p =",g1_v_g2$p.value ,"\n"))
  cat(paste("p<=0.05 vs p>0.01, p =",g1_v_g3$p.value ,"\n"))
}

# ===== Pathogenic promoter and combined ===== 
prom_comb_patho <- bind_rows(gene_muteff$promoter$patho_genes, gene_muteff$combined$patho_genes) %>% distinct()
top_gene_patho_count_prom_comb <- nrow(prom_comb_patho %>% filter(!is.na(cancer)))
top_gene_patho_vec_prom_comb <- c(top_gene_patho_count_prom_comb,  nrow(prom_comb_patho) - top_gene_patho_count_prom_comb)

# use the median
cur_gm_prom_comb <- rbind(
  gene_muteff$promoter$gene_muteff %>% select(gene, min_p_value),
  gene_muteff$combined$gene_muteff %>% select(gene, min_p_value)
) %>% group_by(gene) %>% dplyr::summarise(min_p_value = mean(min_p_value)) %>%
  filter(!(gene%in%prom_comb_patho$gene))

gm_median_prom_com <- median(cur_gm_prom_comb$min_p_value)
gm_ctrl_prom_com <- cur_gm_prom_comb %>% filter(min_p_value > gm_median_prom_com) 
ctrl_gene_patho_count_prom_com <- length(intersect(gm_ctrl_prom_com$gene, pathology$gene))
ctrl_gene_patho_vec_prom_com <- c(ctrl_gene_patho_count_prom_com, length(gm_ctrl_prom_com$gene) - ctrl_gene_patho_count_prom_com)
mat_prom_com <- data.frame(
  "Top genes" = top_gene_patho_vec_prom_comb,
  "Ctrl genes" = ctrl_gene_patho_vec_prom_com,
  row.names = c("Prognostic", "Non-prognostic"),
  stringsAsFactors = FALSE
)
fisher_res <- fisher.test(mat_prom_com, alternative = "greater")

# ==============

# Top gene analyses
top_genes_enhancer <- gene_muteff$enhancer$patho_genes$gene
top_genes_promoter <- gene_muteff$promoter$patho_genes$gene
top_genes_combined <- gene_muteff$combined$patho_genes$gene

tg_promoter_only <- setdiff(top_genes_promoter, top_genes_combined)
length(intersect(pathology$gene, tg_promoter_only))

tg_promoter_only_df <- gene_muteff$promoter$gene_muteff %>% 
    filter(gene %in% tg_promoter_only & p_value_mutated_vs_unmutated < 0.1)

tg_enhancer_df <- gene_muteff$enhancer$gene_muteff %>% 
  filter(gene %in% top_genes_enhancer & p_value_mutated_vs_unmutated < 0.1)


pval_promoter_combined <- (gene_muteff$combined$gene_muteff %>% 
                             mutate("minlogp_combined" = -log(min_adjusted_p_value)) %>%
                             select(c(gene,minlogp_combined,min_adjusted_p_value)) %>%
                             rename(min_adj_p_combined = min_adjusted_p_value)
    ) %>%
      inner_join((
        gene_muteff$promoter$gene_muteff %>% 
          mutate("minlogp_promoter" = -log(min_adjusted_p_value)) %>%
          select(c(gene,minlogp_promoter,min_adjusted_p_value)) %>%
          rename(min_adj_p_promoter = min_adjusted_p_value)
      )) 
  # inner_join((
  #   gene_muteff$enhancer$gene_muteff %>% 
  #     mutate("minlogp_enhancer" = -log(min_p_value)) %>%
  #     select(c(gene,minlogp_enhancer,-min_p_value))
  # ))

write.csv(pval_promoter_combined %>% distinct(), "p_adj_promoter_vs_combined.csv", quote=FALSE, row.names=FALSE)
ggplot(pval_promoter_combined, aes(x=minlogp_combined, y=minlogp_promoter)) +
  geom_point(size=0.5) +
  stat_cor(label.y = 17, 
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_regline_equation(label.y = 16)
