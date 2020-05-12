tf_name <- "MYB"
filename_enhancer <- paste0(tf_name, "_binding_enhancer_p_value.csv")
filename_promoter <- paste0(tf_name, "_binding_promoter_p_value.csv")

# import the p-value result
setwd("/Users/jz132/Desktop/hardac-xfer/clean_data")
data_enhancers_mutated_effect <- read_csv(filename_enhancer)
data_promoters_mutated_effect <- read_csv(filename_promoter)

# annotate the genes associated with each enhancer
setwd("/Users/jz132/r_projects/cancer-mutations/pelinks")
data_eplinks_short <- read_delim("eplinks-fantom-filtered.csv", delim = ",")
data_eplinks_long <- data_eplinks_short %>%
  mutate(tss = strsplit(tss, ";")) %>%
  unnest(tss)

data_enhancers_mutated_effect_annotated <- data_enhancers_mutated_effect %>% 
  inner_join(data_eplinks_short) %>% 
  select(-c("tss", "length", "enhancer")) %>% 
  distinct() %>% 
  mutate(reg_region = "enhancer")

# annotate the genes associated with each promoter
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

data_promoters_mutated_effect_annotated <- data_promoters_mutated_effect %>% 
  mutate(pos = end - 1000) %>%
  inner_join(mapping_refseq_tss) %>%
  inner_join(mapping_hgnc) %>%
  select(-c("pos", "tss", "length", "promoter")) %>%
  distinct() %>% 
  mutate(reg_region = "promoter")

# integrate promoters and enhancers by gene, considering only those elements with mutations
data_reg_mutated_effect_annotated <- bind_rows(data_enhancers_mutated_effect_annotated,
                                               data_promoters_mutated_effect_annotated) %>%
  select(chromosome, start, end, gene, count, p_value, reg_region) %>% 
  arrange(p_value)

data_reg_mutated_effect_summarised <- data_reg_mutated_effect_annotated %>%
  group_by(gene) %>%
  summarise(count = sum(count),
            test_stat = -2*sum(log(p_value)),
            df = 2*n()) %>%
  mutate(p_value = pchisq(test_stat, df = df, lower.tail = F)) %>%
  arrange(p_value)

# refer to cancer gene census
setwd("/Users/jz132/Desktop/cosmic/")

data_cosmic_raw <- read_csv("cosmic_cancer_gene_census.csv")
data_cosmic <- data_cosmic_raw %>% 
  select(gene_symbol = `Gene Symbol`, 
         name = Name, 
         location = `Genome Location`, 
         tier = Tier, 
         somatic = Somatic, 
         germline = Germline,
         tumor_types_somatic = `Tumour Types(Somatic)`, 
         tumor_types_germline = `Tumour Types(Germline)`,
         role_in_cancer = `Role in Cancer`,
         mutation_types = `Mutation Types`)

data_reg_cosmic <- data_reg_mutated_effect_summarised %>% 
  left_join(data_cosmic %>% select(gene = gene_symbol, 
                                   tier,
                                   tumor_types_somatic,
                                   role_in_cancer))
  
# output the result
setwd("/Users/jz132/Desktop/hardac-xfer/output")
filename_effect_output <- paste0(tf_name, "_binding_change_result.csv")
write.csv(data_reg_mutated_effect_summarised, filename_effect_output, row.names = F)

filename_cosmic_output <- paste0(tf_name, "_cosmic_result.csv")
write.csv(data_reg_cosmic, filename_cosmic_output, row.names = F)


# # check with gene expression data
# setwd("/Users/jz132/Downloads/")
# data_expr_raw <- read_tsv("exp_seq.LIRI-JP.tsv")
# 
# data_expr <- data_expr_raw %>%
#   filter(grepl("Cancer", analysis_id)) %>%
#   select(icgc_donor_id, gene_id, normalized_read_count)
# 
# data_expr_select <- data_expr %>% 
#   filter(gene_id == gene_select)
# 
# data_expr_select_case <- data_expr_select %>% 
#   filter(icgc_donor_id %in% data_donor_select$icgc_donor_id) %>% 
#   distinct()
# 
# data_expr_select_control <- data_expr_select %>% 
#   filter(!icgc_donor_id %in% data_donor_select$icgc_donor_id) %>% 
#   distinct()
# 
# boxplot(data_expr_select_case$normalized_read_count,
#         data_expr_select_control$normalized_read_count)
