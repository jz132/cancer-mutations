source("~/r_projects/cancer-mutations/ssm-background/trinucleotide_bg_promoter.R")

promoters_w_mutation <- which(data_promoters_mutated$count > 0)

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

# import the promoter annotations
setwd("/Users/jz132/r_projects/cancer-mutations/pelinks/RefSeq")
mapping_hgnc <- read_delim("refseq_to_hgnc.txt", delim = "\t") %>%
  select(gene = `Approved symbol`,
         refseq_id = `RefSeq IDs`,
         refseq_id_ncbi = `RefSeq(supplied by NCBI)`) %>%
  pivot_longer(cols = refseq_id:refseq_id_ncbi, values_to = "tss") %>%
  select(gene, tss) %>%
  na.omit() %>%
  distinct()

mapping_refseq_tss <- read_delim("refseq_TSS_hg19_170929.bed", delim = '\t',
                                 col_names = F) %>%
  select(c(1,2,4)) %>%
  dplyr::rename(chromosome = "X1", pos = "X2", tss = "X4")

# import the p-value and mutation effect results predicted by 666 PBM data sets
setwd("/Users/jz132/Desktop/hardac-xfer/raw_data/promoter")
p_value_filenames <- Sys.glob("p_value_promoter*.txt")
mut_effect_filenames <- Sys.glob("mut_effect_promoter*.txt")
upbm_filenames <- gsub("p_value_promoter_", "", p_value_filenames)
num_files <- length(p_value_filenames)

p_value_matrix_per_promoter <- NULL
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
  
  data_promoters_mutated_effect <- data_promoters_mutated %>% 
    mutate(p_value = p_value_list) %>% 
    mutate(mut_effect = mut_effect_list) %>%
    filter(count > 0)
  
  data_promoters_mutated_effect_annotated <- data_promoters_mutated_effect %>% 
    mutate(pos = end - 1000) %>%
    inner_join(mapping_refseq_tss) %>%
    inner_join(mapping_hgnc) %>%
    select(-c("pos", "tss", "length", "promoter")) %>%
    distinct() %>%
    mutate(z_score = qnorm(1-p_value/2)*sign(mut_effect)) %>% 
    mutate(reg_region = "promoter")
  
  # Combine p-values using Stouffer's Z-score method
  data_promoters_mutated_effect_summarised <- 
    data_promoters_mutated_effect_annotated %>% 
    group_by(gene) %>%
    summarise(count = sum(count),
              z_score = sum(z_score)/sqrt(n())) %>%
    mutate(p_value = 2*(1-pnorm(abs(z_score)))) %>%
    arrange(gene)
  
  p_value_matrix_per_promoter <- cbind(p_value_matrix_per_promoter,
                                       data_promoters_mutated_effect$p_value)
  p_value_matrix_per_gene <- cbind(p_value_matrix_per_gene, 
                                   data_promoters_mutated_effect_summarised$p_value)
}

rownames(p_value_matrix_per_promoter) <- data_promoters_mutated_effect$promoter
colnames(p_value_matrix_per_promoter) <- gsub("p_value_promoter_", "",
                                              p_value_filenames)
rownames(p_value_matrix_per_gene) <- data_promoters_mutated_effect_summarised$gene
colnames(p_value_matrix_per_gene) <- gsub("p_value_promoter_", "",
                                          p_value_filenames)

adjusted_p_value_matrix_per_gene <- apply(p_value_matrix_per_gene, 2,
                                          p.adjust, method = "hochberg")

# promoter level results
p_value_min_per_promoter <- apply(p_value_matrix_per_promoter, 1, min)

data_promoters_mutated_max_effect_per_promoter <- data_promoters_mutated_effect %>% 
  select(promoter, count) %>%
  mutate(min_p_value = p_value_min_per_promoter) %>% 
  arrange(min_p_value)
hist(log10(data_promoters_mutated_max_effect_per_promoter$min_p_value))
summary(log10(data_promoters_mutated_max_effect_per_promoter$min_p_value))

# gene level results 
p_value_min_per_gene <- apply(p_value_matrix_per_gene, 1, min)
adjusted_p_value_min_per_gene <- apply(adjusted_p_value_matrix_per_gene, 1, min)
filename_which_min_per_gene <- p_value_filenames[apply(p_value_matrix_per_gene, 1, which.min)]
upbm_which_min_per_gene <- gsub("p_value_promoter_", "", filename_which_min_per_gene)

data_promoters_mutated_max_effect_per_gene <- data_promoters_mutated_effect_summarised %>%
  select(gene, count) %>% 
  mutate(min_p_value = p_value_min_per_gene,
         min_adjusted_p_value = adjusted_p_value_min_per_gene,
         upbm = upbm_which_min_per_gene) %>% 
  inner_join(mapping_table_per_upbm) %>%
  arrange(min_p_value)
hist(log10(data_promoters_mutated_max_effect_per_gene$min_p_value))
summary(log10(data_promoters_mutated_max_effect_per_gene$min_p_value))

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
data_promoters_mutated_max_effect_per_gene %>% inner_join(data_cosmic)

# check with expression data
setwd("/Users/jz132/Desktop/exp_seq/")
# data_expr_raw <- read_tsv("exp_seq.LIRI-JP.tsv")
data_expr <- data_expr_raw %>%
  select(icgc_donor_id, gene_id, normalized_read_count,
         analysis_id, submitted_sample_id)

gene_select <- "ATXN3"
data_promoters_select <- data_promoters_mutated_effect_annotated %>% 
  filter(gene == gene_select)

donor_select <- data_icgc_wgs_to_join %>% 
  genome_inner_join(data_promoters_select) %>% 
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
# data_promoters_mutated_effect_summarised <- 
#   data_promoters_mutated_effect_annotated %>% 
#   group_by(gene) %>%
#   summarise(count = sum(count),
#             test_stat = -2*sum(log(p_value)),
#             df = 2*n()) %>%
#   mutate(p_value = pchisq(test_stat, df = df, lower.tail = F)) %>%
#   arrange(p_value)

