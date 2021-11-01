library(tidyverse)
library(data.table)
setwd("/Users/vincentiusmartin/Research/CancerMutation")
source("functions.R")
# ============ Entropy ============ 

entropy_df <- fread_sep("dataset/epdata/enhancer_entropy.csv","enhancer") %>%
  mutate(length = end-start)
ggplot(entropy_df, aes(x=entropy, y=length)) + 
  xlab("Entropy") + 
  ylab("Enhancer length") + 
  geom_point()
ggsave("entropy_vs_enhancer_length.pdf")
  
# promoter: 10.5 to 12.5

z_df <- fread("dataset/combined_predictions/enhancer_muteffect.csv") %>%
  select(chromosome,start,end,absmax_z_score) %>%
  mutate(length = end-start,
         abs_zscore = abs(absmax_z_score))
ggplot(z_df, aes(x=abs_zscore, y=length)) + 
  geom_point(size=0.5)

entropy_z <- entropy_df %>% select(chromosome,start,end,entropy) %>%
  inner_join(z_df %>% select(chromosome, start, end, abs_zscore))
ggplot(entropy_z, aes(x=abs_zscore, y=entropy)) + 
  geom_point(size=0.5)


entropy_df %>% filter(length == min(length))
entropy_df %>% filter(length == max(length))

# ============ PLOT MUTATIONAL BURDEN ============ 
promoter_snp_count <- fread("dataset/epdata/promoter_snp_count.csv") %>%
  mutate(type="promoter")
enhancer_snp_count <- fread("dataset/epdata/enhancer_snp_count.csv") %>%
  mutate(type="enhancer")

median(enhancer_snp_count$length)
ggplot(enhancer_snp_count,aes(x = "",y = enhancer_snp_count$length)) + 
  labs(y="Length") +
  geom_violin(width=0.5, fill = "#F8766D") +
  geom_boxplot(width=0.1)
ggsave("enhancer_length_dist.pdf")

all_snp_count <- as.data.frame(
  rbind(promoter_snp_count, enhancer_snp_count) %>% 
    group_by(type, count) %>% tally()
  ) %>%
  group_by(type, count = replace(count, count > 4, ">4")) %>% 
  summarise(n = sum(n), .groups = 'drop') %>% 
  mutate(count = factor(all_snp_count$count, levels = c(0:4, ">4")))
  
ggplot(data=all_snp_count, aes(x=count, y=n, fill=type)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25) +
  labs(y="Number of regulatory regions", x = "Mutation count") +
  theme_minimal() +
  theme(legend.title=element_blank())
ggsave("cis_snp_count.pdf")

# promoter is generally longer than enhancer

# ============ PLOT TF BINDING CHANGE AS -LOG(P) ============ 

# normal
dfe_normal <- fread("dataset/combined_predictions/enhancer_gene_muteffect.csv") %>%
    mutate(difpos = -log(min_adjusted_p_value),
           ynorm = dnorm(difpos)
           ) %>% select(gene, difpos, ynorm, min_adjusted_p_value)
# rank
dfe <- fread("dataset/combined_predictions/enhancer_gene_muteffect_expanded.csv") 
dfp <- fread("dataset/combined_predictions/promoter_gene_muteffect_expanded.csv") 
dfcmb <- fread("dataset/combined_predictions/combined_gene_muteffect_expanded.csv") 
dfcmb_signif_genes <- dfcmb[dfcmb$min_adjusted_p_value < 0.05]$gene
dfe_signif_genes <- dfe[dfe$min_adjusted_p_value < 0.05]$gene
dfp_signif_genes <- dfp[dfp$min_adjusted_p_value < 0.05]$gene
setdiff(dfcmb_signif_genes, union(dfe_signif_genes,dfp_signif_genes))
length(intersect(dfp_signif_genes, dfcmb_signif_genes))

# library("viridis")
unique_genes <- sort(unique(c(dfcmb_signif_genes,dfe_signif_genes,dfp_signif_genes)))
df_genes_colors <- data.frame(gene=unique_genes, colour=as.character(viridis_pal(option = "A")(length(unique_genes))))
dfplt <- dfcmb %>% left_join(df_genes_colors, on="gene")
#dfplt <- mutate(dfplt, colour=ifelse(is.na(dfplt$colour), "black",as.character(dfplt$colour)))

# geom_text(aes(label=gene),hjust=0, vjust=0)  --> for all genes
# ggplot(dfplt, aes(x = difpos, y = ynorm)) +
#   geom_line() +
#   geom_point(colour=dfplt$colour) + 
#   theme_minimal() +
#   # geom_text(data=subset(dfplt, !is.na(colour)),
#   #           aes(difpos,ynorm,label=gene), vjust=-0.5) +
#   scale_x_continuous("-log(p)") +
#   ylab("Density")

pcut <- 0.05
df_plt <- dfcmb %>% mutate(difpos = -log(min_adjusted_p_value)) %>%
  filter((n_paired_case >= 3) & (n_case_cancer >= 1) ) %>%
  arrange(desc(difpos)) %>%
  mutate(
    category = case_when(
      p_value_mutated_vs_unmutated < pcut & p_value_cancer_vs_normal < pcut ~ "both significant",
      p_value_mutated_vs_unmutated < pcut ~ "mutated_vs_unmutated",
      p_value_cancer_vs_normal < pcut ~ "paired_cancer_vs_normal",
      TRUE ~ "notsignificant"
      ),
    color = case_when(
      p_value_mutated_vs_unmutated < pcut & p_value_cancer_vs_normal < pcut ~ "both significant",
      p_value_mutated_vs_unmutated < pcut ~ "mutated_vs_unmutated",
      p_value_cancer_vs_normal < pcut ~ "paired_cancer_vs_normal",
      TRUE ~ "notsignificant"
    ),
    rank = row_number()
    ) %>%
  dplyr::select(gene, difpos, rank, min_adjusted_p_value, category)
  # filter(rank < 500) %>% 
  # left_join(df_genes_colors, on="gene")
#df_plt <- mutate(df_plt, colour=ifelse(is.na(df_plt$colour), "black",as.character(df_plt$colour)))
#signif_colors <- (df_plt %>% filter(colour != "black"))$colour
# loess_mdl <- loess(df_plt$difpos~df_plt$rank)
# df_plt$difsmooth <- predict(loess_mdl)
# colour=signif_colors, 
# data=subset(df_plt, category != "notsignif")
ggplot(df_plt, aes(x = rank, y = difpos, label=category)) +
  geom_line(color="black", lwd=0.5) +
  geom_point(aes(colour = category), show.legend = TRUE) + 
  theme_minimal() +
  ylab("-log(p)") +
  scale_color_manual(values=c("red","green","gray","blue"))
ggsave("combined_plot.pdf",  width = 10, height = 8, dpi = 150)
# just rank

# ============ ANALYZE PROMOTER AND ENHANCER TOP GENES ============ 

enhancer_annotated <- fread("dataset/annotated/enhancer_annotated.csv")
promoter_annotated <- fread("dataset/annotated/promoter_annotated.csv")
promoter_annotated %>% group_by(gene) %>% count() %>% filter(n > 1)

promoter_genes <- unique(promoter_annotated$gene)
enhancer_genes <- unique(enhancer_annotated$gene)
combined_genes <- union(promoter_genes,enhancer_genes)

res_enhancer <- fread("dataset/combined_predictions/enhancer_gene_muteffect_expanded.csv")
res_promoter <- fread("dataset/combined_predictions/promoter_gene_muteffect_expanded.csv")
res_combined <- fread("dataset/combined_predictions/combined_gene_muteffect_expanded.csv")

pcut <- 0.05
top_promoter <- res_promoter[res_promoter$min_adjusted_p_value < 0.05]
top_promoter_genes <- top_promoter$gene
top_enhancer_genes <- res_enhancer[res_enhancer$min_adjusted_p_value < 0.05]$gene
top_combined <- res_combined[res_combined$min_adjusted_p_value < 0.05]
top_combined_genes <- top_combined$gene

length(top_combined_genes)

length(intersect(top_promoter_genes, top_combined_genes))
length(intersect(top_enhancer_genes, top_combined_genes))
setdiff(top_combined_genes, union(top_promoter_genes,top_enhancer_genes))


# ============ ANALYZE PROMOTER AND ENHANCER GENE EXPRESSION ============ 

top_combined[top_combined$p_value_mutated_vs_unmutated < 0.05]
