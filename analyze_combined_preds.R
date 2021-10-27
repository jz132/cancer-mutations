library(tidyverse)
library(data.table)
library(fuzzyjoin)
rename <- dplyr::rename

setwd("/Users/vincentiusmartin/Research/CancerMutation")

# gene_muteff_path <- "dataset/combined_predictions/combined_gene_muteffect.csv"
# cis_annotated_paths <- c("dataset/annotated/promoter_annotated.csv", "dataset/annotated/enhancer_annotated.csv")
# snps_in_cis_paths <- c("dataset/epdata/snps_in_promoter.csv", "dataset/epdata/snps_in_enhancer.csv")
# cistype <- c("promoter", "enhancer")

# gene_muteff_path <- "dataset/combined_predictions/promoter_gene_muteffect.csv"
# cis_annotated_paths <- c("dataset/annotated/promoter_annotated.csv")
# snps_in_cis_paths <- c("dataset/epdata/snps_in_promoter.csv")
# cistype <- c("promoter")

cosmic_path <- "dataset/input/cosmic_cancer_gene_census.csv"
expr_path <- "dataset/input/icgc_exp_LIRI_JP.tsv"

gene_muteff <- fread(gene_muteff_path)
cis_annotated_combined <- lapply(cis_annotated_paths,fread) %>% bind_rows()
snps_in_cis_combined <- lapply(snps_in_cis_paths,fread) %>% bind_rows()

# ======== Annotate each gene by cosmic tier ======== 
cosmic_df <- fread(cosmic_path, select=c("Gene Symbol","Tier")) %>%
  rename(gene = "Gene Symbol", cosmic_tier = "Tier") %>%
  distinct()
gene_muteff <- gene_muteff %>% 
  left_join(cosmic_df %>% select(gene, cosmic_tier), by="gene")

# ======== Annotate each gene by expression data ======== 
expr_df <- fread(expr_path) %>%
  select(icgc_donor_id,
         gene = gene_id,
         norm_count = normalized_read_count,
         analysis_id) 
genes_w_exp <- unique(expr_df$gene)
gene_muteff <- mutate(gene_muteff, expr_avail = ifelse(gene %in% genes_w_exp, 1, 0))

ngenes <- length(gene_muteff$gene)
div_const <- ngenes %/% 100 # for progress report
expr_pval_df <- lapply(seq(ngenes), function(i){
  if (i %% div_const == 0) {
    writeLines(paste0("Progress ",i,"/",ngenes))
  }
  
  cur_gene <- gene_muteff$gene[i]
  cis_w_gene <- cis_annotated_combined %>%
    filter(gene == cur_gene)
    # inner_join(geneseqs) %>%
    # filter(p_less > 0.9)

  # get donors with mutations in the cis-region that linked to the gene of interest
  donor_select <- cis_w_gene %>%
    genome_inner_join(snps_in_cis_combined, by=c("chromosome", "start", "end")) %>%
    pull(icgc_donor_id) %>%
    unique()
  
  # expression of the gene of interest
  expr_select <- expr_df %>%
    filter(gene == cur_gene)
  
  # Look at expression data with the selected donors
  expr_w_donor <- expr_select %>% 
    filter(icgc_donor_id %in% donor_select) 
  
  expr_cancer_case <- expr_w_donor %>%
    filter(grepl("Cancer", analysis_id)) %>%
    distinct(icgc_donor_id, gene, norm_count)
  
  expr_normal_case <- expr_w_donor %>%
    filter(grepl("Liver", analysis_id)) %>%
    distinct(icgc_donor_id, gene, norm_count)
  
  # paired case
  cancer_vs_normal <- expr_cancer_case %>%
    inner_join(expr_normal_case,
               by = c("icgc_donor_id", "gene"))
  p_cancer_vs_normal <- ifelse(nrow(cancer_vs_normal) == 0, 
                               NA,
                               wilcox.test(cancer_vs_normal$norm_count.x,
                                           cancer_vs_normal$norm_count.y,
                                           paired = T)$p.value)
  
  # Look at expression data in the non selected donor (control)
  expr_wo_donor <- expr_select %>% 
    filter(!icgc_donor_id %in% donor_select) 
  
  expr_cancer_ctrl <- expr_wo_donor %>%
    filter(grepl("Cancer", analysis_id)) %>%
    distinct(icgc_donor_id, gene, norm_count)
  
  p_mutated_vs_unmutated  <- ifelse(nrow(expr_cancer_case) == 0 | nrow(expr_cancer_ctrl)==0, 
                                   NA,
                                   wilcox.test(expr_cancer_case$norm_count,
                                               expr_cancer_ctrl$norm_count, "two.sided")$p.value)
  
  # when plotting is needed
  # (expr_cancer_case$norm_count vs expr_cancer_ctrl$norm_count) 
  # or cancer_vs_normal$norm_count.x,cancer_vs_normal$norm_count.y
  # dat1 <- data.table(group="cancer",value=cancer_vs_normal$norm_count.x)
  # dat2 <- data.table(group="normal",value=cancer_vs_normal$norm_count.y)
  # plot_data <- rbind(dat1,dat2)
  # ggplot(plot_data, aes(x=group, y=value, fill=group)) +
  #   geom_boxplot() +
  #   labs(title=paste0(cur_gene," expression in combined analyses"),
  #        y = "normalized read count", x= NULL) +
  #   theme(legend.position = "none")
    #ggsave(paste0(cur_gene,"_promoter.pdf"))
  
  
  # mutated vs unmutated
  
  dat1 <- data.table(group="mutated",value=expr_cancer_case$norm_count)
  dat2 <- data.table(group="unmutated",value=expr_cancer_ctrl$norm_count)
  plot_data <- rbind(dat1,dat2)
  ggplot(plot_data, aes(x=group, y=value, fill=group)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#FFA500","#989898")) +
  labs(title=paste0(cur_gene," expression in combined"),
         y = "normalized read count", x= NULL) +
       ylim(0,47) +
       theme(legend.position = "none")
  ggsave(paste0(cur_gene,"_combined_top_vs_ctrl.pdf"))
  
  
  return (list(p_value_cancer_vs_normal = p_cancer_vs_normal, 
               n_paired_case = nrow(cancer_vs_normal),
               p_value_mutated_vs_unmutated = p_mutated_vs_unmutated,
               n_case_cancer = nrow(expr_cancer_case),
               n_control_cancer = nrow(expr_cancer_ctrl)
               ))
}) %>% bind_rows()

gene_muteff <-  bind_cols(gene_muteff, expr_pval_df)
cisname <- ifelse(length(cistype) == 1, cistype[1], "combined")
write.csv(gene_muteff %>% distinct(), paste0(cisname,"_gene_muteffect_expanded.csv"), quote=FALSE, row.names=FALSE)
