library(tidyverse)
library(data.table)
library(ggplot2)

setwd("/Users/vincentiusmartin/Research/CancerMutation")
preg_df <- fread("p_act_promoter_vs_combined.csv") %>%
  mutate(
    significance = case_when(
      min_adj_p_combined < 0.05 & min_adj_p_promoter < 0.05 ~ "both",
      min_adj_p_combined < 0.05 ~ "combined_only",
      min_adj_p_promoter < 0.05 ~ "promoter_only",
      TRUE                      ~ "notsignif"
    )
  )

# get gene annotated regulatory regions
cis_annotated_promoter <- fread("dataset/annotated/promoter_annotated.csv")
cis_annotated_promoter$cisloc <- "promoter"
cis_annotated_enhancer <- fread("dataset/annotated/enhancer_annotated.csv")
cis_annotated_enhancer$cisloc <- "enhancer"
cis_annotated <- rbind(cis_annotated_promoter,cis_annotated_enhancer) %>% 
    select(gene, cisloc) %>%
    group_by(gene) %>%
    mutate(cisloc = paste0(unique(sort(cisloc)), collapse = "-")) %>%
    distinct()

gene_locs_p <- preg_df %>% 
  inner_join(cis_annotated, by="gene")
promoter_cutoff <- tail((
  gene_locs_p %>% filter(min_adj_p_promoter < 0.05) %>% arrange(min_adj_p_promoter))$minlogp_promoter, n = 1)
combined_cutoff <- tail((
  gene_locs_p %>% filter(min_adj_p_combined < 0.05) %>% arrange(min_adj_p_combined))$minlogp_combined, n = 1)
  
ggplot(gene_locs_p, aes(x=minlogp_combined, y=minlogp_promoter)) +
  geom_point(size=0.5, aes(color=significance)) +
  geom_hline(yintercept=promoter_cutoff) +
  geom_vline(xintercept=combined_cutoff) +
  guides(color = guide_legend(override.aes = list(size = 3) ) )
  #geom_text(aes(label=ifelse(min_adj_p_combined<0.05|min_adj_p_promoter<0.05,as.character(gene),'')),size=2,hjust=0,vjust=0)+
# stat_cor(label.y = 18, 
#            aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
#   stat_regline_equation(label.y = 17) 
ggsave("scatter_p_act_promoter_vs_combined.pdf")
      