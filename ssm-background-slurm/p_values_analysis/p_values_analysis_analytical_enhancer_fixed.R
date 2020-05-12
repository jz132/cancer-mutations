source("~/r_projects/cancer-mutations/ssm-background/trinucleotide_bg_enhancer.R")

enhancers_w_mutation <- which(data_enhancers_mutated$count > 0)

# import the mapping between TFs and PBMs
# there are two naming issues I am still discussing with Vincentius
setwd("/Users/jz132/Desktop/hardac-xfer/mappings")
mapping_table_per_tf <- read_csv("mapping_to_12mer_tbl.csv") %>% 
  select(tf = gene, 
         upbm) %>% 
  mutate(upbm = gsub(" ", "", upbm, fixed = T)) %>% 
  mutate(upbm = gsub("prediction6mer.Mus_musculus|M01338_1.94d|Badis09|Six6_2267.4=v1.txt",
                     "prediction6mer.Mus_musculus|M01339_1.94d|Berger08|Six6_2267.4.txt",
                     upbm, fixed = T))

mapping_table_per_upbm <- mapping_table_per_tf %>% 
  mutate(upbm = strsplit(upbm, ";")) %>% 
  unnest(upbm) %>% 
  group_by(upbm) %>% 
  summarise(tf = paste0(tf, collapse = ";")) %>%
  filter(upbm != "prediction6mer.Homo_sapiens|NA|MartinZhao2019|RELA.txt")

# import the enhancer annotations
setwd("/Users/jz132/r_projects/cancer-mutations/pelinks")
data_eplinks_short <- read_delim("eplinks-fantom-filtered.csv", delim = ",")
data_eplinks_long <- data_eplinks_short %>%
  mutate(tss = strsplit(tss, ";")) %>%
  unnest(tss)

# import the p-value and mutation effect results predicted by 666 PBM data sets
setwd("/Users/jz132/Desktop/hardac-xfer/raw_data/enhancer")
p_value_filenames <- Sys.glob("p_value_enhancer*.txt")
mut_effect_filenames <- Sys.glob("mut_effect_enhancer*.txt")
upbm_filenames <- gsub("p_value_enhancer_", "", p_value_filenames)
num_files <- length(p_value_filenames)

p_value_matrix_per_enhancer <- NULL
p_value_matrix_per_gene <- NULL

for(i in 1:num_files){
  p_value_filename <- p_value_filenames[i]
  p_value_list <- read_delim(p_value_filename, delim = " ", 
                             col_names = "p_value") %>% 
    pull(p_value)
  mut_effect_filename <- mut_effect_filenames[i]
  mut_effect_list <- read_delim(mut_effect_filename, delim = " ", 
                                col_names = "mut_effect") %>% 
    pull(mut_effect)

  data_enhancers_mutated_effect <- data_enhancers_mutated %>% 
    mutate(p_value = p_value_list) %>% 
    mutate(mut_effect = mut_effect_list) %>%
    filter(count > 0)

  data_enhancers_mutated_effect_annotated <- data_enhancers_mutated_effect %>% 
    inner_join(data_eplinks_short) %>% 
    select(-c("tss", "length", "enhancer")) %>% 
    distinct() %>% 
    mutate(z_score = qnorm(1-p_value/2)*sign(mut_effect)) %>% 
    mutate(reg_region = "enhancer")
  
  # Combine p-values using Stouffer's Z-score method
  data_enhancers_mutated_effect_summarised <- 
    data_enhancers_mutated_effect_annotated %>% 
    group_by(gene) %>%
    summarise(count = sum(count),
              z_score = sum(z_score)/sqrt(n())) %>%
    mutate(p_value = 2*(1-pnorm(abs(z_score)))) %>%
    arrange(gene)
  
  p_value_matrix_per_enhancer <- cbind(p_value_matrix_per_enhancer,
                                       data_enhancers_mutated_effect$p_value)
  p_value_matrix_per_gene <- cbind(p_value_matrix_per_gene, 
                          data_enhancers_mutated_effect_summarised$p_value)
}

rownames(p_value_matrix_per_enhancer) <- data_enhancers_mutated_effect$enhancer
colnames(p_value_matrix_per_enhancer) <- gsub("p_value_enhancer_", "",
                                              p_value_filenames)
rownames(p_value_matrix_per_gene) <- data_enhancers_mutated_effect_summarised$gene
colnames(p_value_matrix_per_gene) <- gsub("p_value_enhancer_", "",
                                          p_value_filenames)

adjusted_p_value_matrix_per_gene <- apply(p_value_matrix_per_gene, 2,
                                          p.adjust, method = "hochberg")

# enhancer level results
p_value_min_per_enhancer <- apply(p_value_matrix_per_enhancer, 1, min)

data_enhancers_mutated_max_effect_per_enhancer <- data_enhancers_mutated_effect %>% 
  select(enhancer, count) %>%
  mutate(min_p_value = p_value_min_per_enhancer) %>% 
  arrange(min_p_value)
hist(log10(data_enhancers_mutated_max_effect_per_enhancer$min_p_value))
summary(log10(data_enhancers_mutated_max_effect_per_enhancer$min_p_value))

# gene level results 
p_value_min_per_gene <- apply(p_value_matrix_per_gene, 1, min)
adjusted_p_value_min_per_gene <- apply(adjusted_p_value_matrix_per_gene, 1, min)
filename_which_min_per_gene <- p_value_filenames[apply(p_value_matrix_per_gene, 1, which.min)]
upbm_which_min_per_gene <- gsub("p_value_enhancer_", "", filename_which_min_per_gene)

data_enhancers_mutated_max_effect_per_gene <- data_enhancers_mutated_effect_summarised %>%
  select(gene, count) %>% 
  mutate(min_p_value = p_value_min_per_gene,
         min_adjusted_p_value = adjusted_p_value_min_per_gene,
         upbm = upbm_which_min_per_gene) %>% 
  inner_join(mapping_table_per_upbm) %>%
  arrange(min_p_value)
hist(log10(data_enhancers_mutated_max_effect_per_gene$min_p_value))
summary(log10(data_enhancers_mutated_max_effect_per_gene$min_p_value))

# check with COSMIC genes
setwd("/Users/jz132/Desktop/cosmic/")
data_cosmic_raw <- read_csv("cosmic_cancer_gene_census.csv")
data_cosmic <- data_cosmic_raw %>% 
  select(gene = `Gene Symbol`, 
         name = Name, 
         location = `Genome Location`, 
         tier = Tier, 
         somatic = Somatic, 
         germline = Germline,
         tumor_types_somatic = `Tumour Types(Somatic)`, 
         tumor_types_germline = `Tumour Types(Germline)`,
         role_in_cancer = `Role in Cancer`,
         mutation_types = `Mutation Types`) 
data_enhancers_mutated_max_effect_per_gene %>% inner_join(data_cosmic)

# check with expression data
setwd("/Users/jz132/Desktop/exp_seq/")
# data_expr_raw <- read_tsv("exp_seq.LIRI-JP.tsv")
data_expr <- data_expr_raw %>%
  select(icgc_donor_id, gene_id, normalized_read_count,
         analysis_id, submitted_sample_id)

gene_select <- "ATXN3"
data_enhancers_select <- data_enhancers_mutated_effect_annotated %>% 
  filter(gene == gene_select)

donor_select <- data_icgc_wgs_to_join %>% 
  genome_inner_join(data_enhancers_select) %>% 
  pull(icgc_donor_id)

data_expr_select <- data_expr %>% 
  filter(gene_id == gene_select)

data_expr_select_case_cancer <- data_expr_select %>%
  filter(icgc_donor_id %in% donor_select) %>%
  filter(grepl("Cancer", analysis_id)) %>%
  select(icgc_donor_id, gene_id, normalized_read_count) %>% 
  distinct()

data_expr_select_case_liver <- data_expr_select %>%
  filter(icgc_donor_id %in% donor_select) %>%
  filter(grepl("Liver", analysis_id)) %>%
  select(icgc_donor_id, gene_id, normalized_read_count) %>% 
  distinct()

data_expr_select_control_cancer <- data_expr_select %>%
  filter(!icgc_donor_id %in% donor_select) %>% 
  filter(grepl("Cancer", analysis_id)) %>%
  select(icgc_donor_id, gene_id, normalized_read_count) %>% 
  distinct()

data_expr_select_control_liver <- data_expr_select %>%
  filter(!icgc_donor_id %in% donor_select) %>% 
  filter(grepl("Liver", analysis_id)) %>%
  select(icgc_donor_id, gene_id, normalized_read_count) %>% 
  distinct()

boxplot(data_expr_select_case_cancer$normalized_read_count,
        data_expr_select_case_liver$normalized_read_count,
        data_expr_select_control_cancer$normalized_read_count,
        data_expr_select_control_liver$normalized_read_count,
        main = paste(gene_select, "expression"),
        names = c("case_cancer", "case_liver", "control_cancer", "control_liver"))

# # Using Fisher's method is not appropriate here as the sign of mut_effect matters
# data_enhancers_mutated_effect_summarised <- 
#   data_enhancers_mutated_effect_annotated %>% 
#   group_by(gene) %>%
#   summarise(count = sum(count),
#             test_stat = -2*sum(log(p_value)),
#             df = 2*n()) %>%
#   mutate(p_value = pchisq(test_stat, df = df, lower.tail = F)) %>%
#   arrange(p_value)


# # additional analysis on whether different tables predict similarly for corresponding tf
# tf_names <- mapping_table$gene
# for(tf_name in tf_names){
#   upbm <- mapping_table %>%
#     filter(gene == tf_name) %>%
#     pull(upbm)
#   
#   upbm_filenames <- unlist(strsplit(upbm, split = ";"))
#   p_value_filenames <- paste0("p_value_enhancer_", upbm_filenames)
#   mut_effect_filenames <- paste0("mut_effect_enhancer_", upbm_filenames)
#   num_files <- length(p_value_filenames)
#   
#   data_enhancers_mutated_effect_list <- list()
#   setwd("/Users/jz132/Desktop/hardac-xfer/raw_data/enhancer")
#   for(i in 1:num_files){
#     p_value_filename <- p_value_filenames[i]
#     p_value_list <- read_delim(p_value_filename, delim = " ", 
#                                col_names = "p_value") %>% 
#       pull(p_value)
#     
#     mut_effect_filename <- mut_effect_filenames[i]
#     mut_effect_list <- read_delim(mut_effect_filename, delim = " ", 
#                                   col_names = "mut_effect") %>% 
#       pull(mut_effect)
#     
#     data_enhancers_mutated_effect <- data_enhancers_mutated %>% 
#       mutate(p_value = p_value_list) %>% 
#       mutate(mut_effect = mut_effect_list) %>%
#       filter(count > 0)
#     
#     data_enhancers_mutated_effect_list[[i]] <- data_enhancers_mutated_effect
#   }
#   
#   names(data_enhancers_mutated_effect_list) <- upbm_filenames
#   print(data_enhancers_mutated_effect_list)
#   
#   cor_plot_names <- unlist(lapply(strsplit(upbm_filenames, split = "|", fixed = T), "[", 4))
#   
#   if(length(data_enhancers_mutated_effect_list) > 1){
#     table_mut_effect <- matrix(0, nrow = length(enhancers_w_mutation),
#                                ncol = length(data_enhancers_mutated_effect_list))
#     for(j in 1:length(data_enhancers_mutated_effect_list)){
#       table_mut_effect[,j] <- data_enhancers_mutated_effect_list[[j]]$mut_effect
#     }
#     colnames(table_mut_effect) <- cor_plot_names
#     setwd("/Users/jz132/Desktop/hardac-xfer/plots")
#     png(file = paste0(tf_name, ".png"), width = 1600, height = 1600, res = 160)
#     pairs.panels(table_mut_effect, digits = 4, smooth = F, ellipses = F,
#                  main = tf_name)
#     dev.off()
#   }
# }
