library(Biostrings)
library(GenomicRanges)
library(tidyverse)
library(fuzzyjoin)

# choose the ICGC study we want to analyze
datasetname <- "LIRI-JP"

# file paths and file names
icgc.data.path <- "~/Desktop/Gordanlab/Data/ICGC"
genomic.interval.path <- "~/r_projects/cancer-mutations/pelinks"
entropy.path <- "~/r_projects/cancer-mutations/pelinks/entropy"
tf.pbm.mapping.path <- "~/Desktop/Gordanlab/Data/qbic/mapping"
jpred.result.path <- paste0("~/Desktop/Gordanlab/Data/jpred/", datasetname)
cosmic.path <- "~/Desktop/Gordanlab/Data/cosmic"
expression.path <- "~/Desktop/Gordanlab/Data/ICGC/exp_seq"
output.path <- "~/r_projects/cancer-mutations/ssm-background-slurm/output"

filename_icgc_ssm <- paste0("simple_somatic_mutation.open.", datasetname, ".tsv")
filename_icgc_exp <- paste0("exp_seq.", datasetname, ".tsv")
filename_enhancer_entropy <- paste0("enhancers_entropy_", datasetname, ".csv")
filename_promoter_entropy <- paste0("promoters_entropy_", datasetname, ".csv")

keep_cols <- c("icgc_mutation_id", 
               "icgc_donor_id", 
               "chromosome",
               "chromosome_start",
               "chromosome_end",
               "assembly_version",
               "mutation_type",
               "reference_genome_allele",
               "mutated_from_allele",
               "mutated_to_allele",
               "consequence_type",
               "gene_affected",
               "transcript_affected",
               "sequencing_strategy")

# import genomic coordinates of promoter/enhancer and links between them and tss
setwd(genomic.interval.path)
data_exon_refseq <- read_delim("all_exons_refseq.txt", delim = "\t", 
                                col_names = c("chromosome", "start", "end", "exon"))

data_promoter_refseq <- read_delim("all_promoters_refseq.txt", delim = "\t",
                                    col_names = c("chromosome", "start", "end", "promoter")) %>%
  arrange(chromosome, start, end)

data_enhancer_fantom <- read_delim("all_enhancers_fantom.txt", delim = "\t",
                                    col_names = c("chromosome", "start", "end", "enhancer")) %>%
  arrange(chromosome, start, end)

mapping_hgnc <- read_delim("refseq_to_hgnc.txt", delim = "\t") %>%
  select(gene = `Approved symbol`,
         refseq_id = `RefSeq IDs`,
         refseq_id_ncbi = `RefSeq(supplied by NCBI)`) %>%
  pivot_longer(cols = refseq_id:refseq_id_ncbi, values_to = "tss") %>%
  select(gene, tss) %>%
  filter(grepl("NM", tss)) %>%
  na.omit() %>%
  distinct()

mapping_ensembl <- read_delim("ensembl_to_hgnc.txt", delim = "\t") %>% 
  select(gene = `Approved symbol`,
         ensembl_id = `Ensembl ID(supplied by Ensembl)`) %>% 
  na.omit() %>% 
  distinct()

mapping_refseq_tss <- read_delim("refseq_TSS_hg19_170929.bed", delim = '\t',
                                 col_names = F) %>%
  select(c(1,2,4)) %>%
  dplyr::rename(chromosome = "X1", pos = "X2", tss = "X4") %>% 
  filter(grepl("NM", tss))

data_eplinks_short <- read_delim("eplinks-fantom-filtered.csv", delim = ",")
data_eplinks_long <- data_eplinks_short %>%
  mutate(tss = strsplit(tss, ";")) %>%
  unnest(tss)

# import enhancer and promoter entropy for weighting
setwd(entropy.path)
data_enhancer_entropy <- read_csv(filename_enhancer_entropy)
data_promoter_entropy <- read_csv(filename_promoter_entropy)

# import ICGC data
setwd(icgc.data.path)
data_icgc_try <- read_tsv(filename_icgc_ssm, col_names = T, n_max = 5)
col_indicator <- paste(ifelse(colnames(data_icgc_try) %in% keep_cols, "?", "-"), collapse = "")
data_icgc_raw <- read_tsv(filename_icgc_ssm, col_names = T, col_types = col_indicator) %>%
  mutate(chromosome = paste0("chr", chromosome))

# focus on single base substitution in WGS only
data_icgc_wgs <- data_icgc_raw %>% 
  filter(sequencing_strategy == "WGS", mutation_type == "single base substitution")
num_donors <- length(data_icgc_wgs %>% distinct(icgc_donor_id) %>% pull())
num_mutations <- length(data_icgc_wgs %>% distinct(icgc_mutation_id, icgc_donor_id) %>% pull())
wgs_donors <- unique(data_icgc_wgs$icgc_donor_id)
if(!all(data_icgc_wgs$mutated_from_allele == data_icgc_wgs$reference_genome_allele)){
  warning("reference genome allele is not the same as mutated from allele")
}

# keep regulatory mutations only
# filter out mutations that are not in enhancer
data_mutations_enhancer <- data_icgc_wgs %>%
  filter(!is.na(consequence_type)) %>% 
  select(chromosome = chromosome,
         start = chromosome_start,
         end = chromosome_end,
         icgc_mutation_id,
         icgc_donor_id,
         ref = reference_genome_allele,
         mut = mutated_to_allele) %>%
  genome_anti_join(data_exon_refseq, by = c("chromosome", "start", "end")) %>%
  genome_inner_join(data_enhancer_fantom, by = c("chromosome", "start", "end")) %>%
  select(chromosome = chromosome.x,
         start = start.x,
         end = end.x,
         icgc_mutation_id, icgc_donor_id, ref, mut) %>%
  distinct()

# filter out mutations that are not in promoter
data_mutations_promoter <- data_icgc_wgs %>%
  filter(!is.na(consequence_type)) %>% 
  select(chromosome = chromosome,
         start = chromosome_start,
         end = chromosome_end,
         icgc_mutation_id,
         icgc_donor_id,
         ref = reference_genome_allele,
         mut = mutated_to_allele) %>%
  genome_anti_join(data_exon_refseq, by = c("chromosome", "start", "end")) %>%
  genome_inner_join(data_promoter_refseq, by = c("chromosome", "start", "end")) %>%
  select(chromosome = chromosome.x,
         start = start.x,
         end = end.x,
         icgc_mutation_id, icgc_donor_id, ref, mut) %>%
  distinct()

data_mutations_reg <- bind_rows(data_mutations_promoter,
                                data_mutations_enhancer)

# for each enhancer, count the number of mutations in it
data_enhancer_mutated <- data_enhancer_fantom %>% 
  genome_left_join(data_mutations_enhancer, by = c("chromosome", "start", "end")) %>%
  select(chromosome = chromosome.x,
         start = start.x,
         end = end.x,
         enhancer,
         icgc_mutation_id,
         icgc_donor_id) %>%
  mutate(count = ifelse(is.na(icgc_mutation_id), 0, 1)) %>%
  group_by(chromosome, start, end, enhancer) %>%
  summarise(count = sum(count)) %>%
  mutate(length = end - start + 1) %>%
  ungroup() %>% 
  inner_join(data_enhancer_entropy)

# for each promoter, count the number of mutations in it
data_promoter_mutated <- data_promoter_refseq %>% 
  genome_left_join(data_mutations_promoter, by = c("chromosome", "start", "end")) %>%
  select(chromosome = chromosome.x,
         start = start.x,
         end = end.x,
         promoter,
         icgc_mutation_id,
         icgc_donor_id) %>%
  mutate(count = ifelse(is.na(icgc_mutation_id), 0, 1)) %>%
  group_by(chromosome, start, end, promoter) %>%
  summarise(count = sum(count)) %>%
  mutate(length = end - start + 1) %>%
  ungroup() %>% 
  inner_join(data_promoter_entropy)

# annotate the genes each promoter is associated with
# the inner join with mapping_hgnc will lead to loss of promoter without mapped gene
data_enhancer_annotated <- data_enhancer_mutated %>%
  inner_join(data_eplinks_short) %>% 
  select(-tss) %>% 
  distinct()

# annotate the genes each enhancer is associated with
# enhancer not linked to any gene are excluded
data_promoter_annotated <- data_promoter_mutated %>%
  mutate(pos = end - 1000) %>%
  inner_join(mapping_refseq_tss) %>%
  inner_join(mapping_hgnc) %>%
  select(-c("pos", "tss")) %>%
  distinct()

data_reg_annotated <- bind_rows(data_promoter_annotated %>% select(-promoter),
                                data_enhancer_annotated %>% select(-enhancer))

# import the mapping between TFs and PBMs
# there are two naming issues I need to discuss with Vincentius
setwd(tf.pbm.mapping.path)
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

# import the p-value and mutation effect results predicted by 667 PBM data sets
setwd(jpred.result.path)
jpred_enhancer_filenames <- Sys.glob("enhancer*.txt")
jpred_promoter_filenames <- Sys.glob("promoter*.txt")
upbm_filenames <- gsub(".*prediction", "prediction", jpred_enhancer_filenames)
num_files <- length(jpred_enhancer_filenames)

z_score_matrix_per_gene <- NULL

for(i in 1:num_files){
  # promoter part
  jpred_promoter_filename <- jpred_promoter_filenames[i]
  
  data_jpred_promoter <- read_delim(jpred_promoter_filename, delim = " ")
  
  data_promoter_mutated_effect <- data_promoter_mutated %>% 
    inner_join(data_jpred_promoter) %>% 
    mutate(z_score = qnorm(p_less)) %>% 
    mutate(p_value = 2*(1-pnorm(abs(z_score))))
  
  data_promoter_mutated_effect_annotated <- data_promoter_mutated_effect %>% 
    mutate(pos = end - 1000) %>%
    inner_join(mapping_refseq_tss) %>%
    inner_join(mapping_hgnc) %>%
    mutate(reg_region = "promoter") %>%
    select(-c("pos", "tss", "length", "promoter", "p_value")) %>% 
    distinct()
  
  # enhancer part
  jpred_enhancer_filename <- jpred_enhancer_filenames[i]
  
  data_jpred_enhancer <- read_delim(jpred_enhancer_filename, delim = " ")
  
  data_enhancer_mutated_effect <- data_enhancer_mutated %>% 
    inner_join(data_jpred_enhancer) %>% 
    mutate(z_score = qnorm(p_less)) %>% 
    mutate(p_value = 2*(1-pnorm(abs(z_score))))
  
  data_enhancer_mutated_effect_annotated <- data_enhancer_mutated_effect %>% 
    inner_join(data_eplinks_short) %>% 
    mutate(reg_region = "enhancer") %>%
    select(-c("tss", "length", "enhancer", "p_value")) %>% 
    distinct()
  
  # integrate promoter and enhancer by gene
  data_reg_mutated_effect_annotated <- bind_rows(data_enhancer_mutated_effect_annotated,
                                                 data_promoter_mutated_effect_annotated) %>%
    select(chromosome, start, end, gene, count, p_less, z_score, s_entropy, reg_region) 
  
  # Combine p-values using Stouffer's Z-score method
  data_reg_mutated_effect_summarised <- 
    data_reg_mutated_effect_annotated %>% 
    group_by(gene) %>%
    summarise(count = sum(count),
              z_score = sum(s_entropy*z_score)/sqrt(sum(s_entropy^2))) %>%
    mutate(p_value = 2*(1-pnorm(abs(z_score)))) %>%
    arrange(gene)
  
  z_score_matrix_per_gene <- cbind(z_score_matrix_per_gene, 
                                   data_reg_mutated_effect_summarised$z_score)
}

rownames(z_score_matrix_per_gene) <- data_reg_mutated_effect_summarised$gene
colnames(z_score_matrix_per_gene) <- upbm_filenames

p_value_matrix_per_gene <- 2*(1-pnorm(abs(z_score_matrix_per_gene)))

# adjusted p value after correcting for multiple comparison
adjusted_p_value_matrix_per_gene <- apply(p_value_matrix_per_gene, 2,
                                          p.adjust, method = "hochberg")

# gene level results 
absmax <- function(x) {x[which.max(abs(x))]}
z_score_absmax_per_gene<- apply(z_score_matrix_per_gene, 1, absmax)
p_value_min_per_gene <- apply(p_value_matrix_per_gene, 1, min)
adjusted_p_value_min_per_gene <- apply(adjusted_p_value_matrix_per_gene, 1, min)
upbm_which_min_per_gene <- upbm_filenames[apply(p_value_matrix_per_gene, 1, 
                                                which.min)]

data_reg_mutated_max_effect_per_gene <- data_reg_mutated_effect_summarised %>%
  select(gene, count) %>% 
  mutate(absmax_z_score = z_score_absmax_per_gene,
         min_p_value = p_value_min_per_gene,
         min_adjusted_p_value = adjusted_p_value_min_per_gene,
         upbm = upbm_which_min_per_gene) %>% 
  inner_join(mapping_table_per_upbm) %>%
  arrange(min_adjusted_p_value)
hist(log10(data_reg_mutated_max_effect_per_gene$min_p_value))
summary(log10(data_reg_mutated_max_effect_per_gene$min_p_value))

# annotate each gene by COSMIC tier
setwd(cosmic.path)
data_cosmic_raw <- read_csv("cosmic_cancer_gene_census.csv")
data_cosmic <- data_cosmic_raw %>% 
  select(gene = `Gene Symbol`, 
         name = Name, 
         location = `Genome Location`, 
         cosmic_tier = Tier, 
         somatic = Somatic, 
         germline = Germline,
         tumor_types_somatic = `Tumour Types(Somatic)`, 
         tumor_types_germline = `Tumour Types(Germline)`,
         role_in_cancer = `Role in Cancer`,
         mutation_types = `Mutation Types`) 
data_reg_result_per_gene <- data_reg_mutated_max_effect_per_gene %>% 
  left_join(data_cosmic %>% select(gene, cosmic_tier))

setwd(output.path)
filename_output <- paste0(datasetname, "_reg_result_brief.csv")
write.csv(data_reg_result_per_gene, filename_output, row.names = F, quote = F)

# annotate each gene by expression data
setwd(expression.path)
if(length(Sys.glob(filename_icgc_exp)) == 0){
  data_reg_result_per_gene <- data_reg_result_per_gene %>%
    mutate(expr_avail = 0)
} else{
  data_expr_raw <- read_tsv(filename_icgc_exp)
  data_expr <- data_expr_raw %>%
    select(icgc_donor_id,
           gene = gene_id,
           norm_count = normalized_read_count,
           analysis_id,
           submitted_sample_id) %>% 
    filter(icgc_donor_id %in% wgs_donors)
  
  if(all(substring(data_expr$gene, 1, 4) == "ENSG")){
    data_expr <- data_expr %>%
      mutate(gene = sapply(strsplit(gene, ".", fixed = T), "[[", 1)) %>%
      rename(ensembl_id = gene) %>%
      inner_join(mapping_ensembl) %>%
      select(-ensembl_id)
  }
  
  genes_w_expression <- unique(data_expr$gene)
  
  data_reg_result_per_gene <- 
    data_reg_result_per_gene %>%
    mutate(expr_avail = ifelse(gene %in% genes_w_expression, 1, 0))
  
  # for donors w mutations, test their cancer cell against their normal cell
  test_cancer_vs_normal <- NULL
  vec_n_paired_case <- NULL
  # for cancer cell, test donors w mutations against donors without mutations
  test_mutated_vs_unmutated <- NULL
  vec_n_case_cancer <- NULL
  vec_n_control_cancer <- NULL
  
  for(gene_select in data_reg_result_per_gene$gene){
    data_reg_annotated_select <- data_reg_annotated %>%
      filter(gene == gene_select)
    
    donor_select <- data_mutations_reg %>%
      genome_inner_join(data_reg_annotated_select) %>%
      pull(icgc_donor_id) %>%
      unique()
    
    data_expr_select <- data_expr %>%
      filter(gene == gene_select)
    
    data_expr_select_case_cancer <- data_expr_select %>%
      filter(icgc_donor_id %in% donor_select) %>%
      filter(grepl("Cancer", analysis_id)) %>%
      select(icgc_donor_id, gene, norm_count) %>%
      distinct()
    
    data_expr_select_case_normal <- data_expr_select %>%
      filter(icgc_donor_id %in% donor_select) %>%
      filter(grepl("Liver", analysis_id)) %>%
      select(icgc_donor_id, gene, norm_count) %>%
      distinct()
    
    data_expr_select_control_cancer <- data_expr_select %>%
      filter(!icgc_donor_id %in% donor_select) %>%
      filter(grepl("Cancer", analysis_id)) %>%
      select(icgc_donor_id, gene, norm_count) %>%
      distinct()
    
    data_expr_select_control_normal <- data_expr_select %>%
      filter(!icgc_donor_id %in% donor_select) %>%
      filter(grepl("Liver", analysis_id)) %>%
      select(icgc_donor_id, gene, norm_count) %>%
      distinct()
    
    data_test_cancer_vs_normal <- data_expr_select_case_cancer %>%
      inner_join(data_expr_select_case_normal,
                 by = c("icgc_donor_id", "gene"))
    n_paired_case <- nrow(data_test_cancer_vs_normal)
    
    if(n_paired_case==0){
      p_cancer_vs_normal <- NA
    } else{
      p_cancer_vs_normal <- wilcox.test(data_test_cancer_vs_normal$norm_count.x,
                                        data_test_cancer_vs_normal$norm_count.y,
                                        paired = T
      )$p.value
    }
    
    n_case_cancer <- nrow(data_expr_select_case_cancer)
    n_control_cancer <- nrow(data_expr_select_control_cancer)
    
    if(n_case_cancer==0||n_control_cancer==0){
      p_mutated_vs_unmutated <- NA
    } else{
      p_mutated_vs_unmutated <- wilcox.test(data_expr_select_case_cancer$norm_count,
                                            data_expr_select_control_cancer$norm_count
      )$p.value
    }
    
    test_cancer_vs_normal <- c(test_cancer_vs_normal, p_cancer_vs_normal)
    vec_n_paired_case <- c(vec_n_paired_case, n_paired_case)
    test_mutated_vs_unmutated <- c(test_mutated_vs_unmutated, p_mutated_vs_unmutated)
    vec_n_case_cancer <- c(vec_n_case_cancer, n_case_cancer)
    vec_n_control_cancer <- c(vec_n_control_cancer, n_control_cancer)
  }
  
  data_reg_result_per_gene <- data_reg_result_per_gene %>%
    mutate(p_value_cancer_vs_normal = test_cancer_vs_normal,
           n_paired_case = vec_n_paired_case,
           p_value_mutated_vs_unmutated = test_mutated_vs_unmutated,
           n_case_cancer = vec_n_case_cancer,
           n_control_cancer = vec_n_control_cancer)
}

setwd(output.path)
filename_output <- paste0(datasetname, "_reg_result_per_gene_weighted.csv")
write.csv(data_reg_result_per_gene, filename_output, row.names = F, quote = F)

# expand the result to include all significant TFs for each top gene
p_cutoff <- 0.1
data_reg_result_per_gene_top <- data_reg_result_per_gene %>%
  filter(min_adjusted_p_value < p_cutoff)

top_genes <- data_reg_result_per_gene_top %>% pull(gene)

data_reg_result_per_gene_top_expanded <- NULL
for(gene in top_genes){
  p_values <- p_value_matrix_per_gene[gene,]
  adjusted_p_values <- adjusted_p_value_matrix_per_gene[gene,]
  
  upbm_select <- names(adjusted_p_values[adjusted_p_values < p_cutoff])
  adjusted_p_values_select <- adjusted_p_values[upbm_select]
  p_values_select <- p_values[upbm_select]
  
  temp_result <- tibble(gene = gene,
                        p_value = p_values_select,
                        adjusted_p_value = adjusted_p_values_select,
                        upbm = upbm_select) %>%
    inner_join(data_reg_result_per_gene_top %>%
                 select(-c("min_p_value", "min_adjusted_p_value",
                           "upbm", "tf"))) %>%
    inner_join(mapping_table_per_upbm) %>%
    select(gene, count, p_value, adjusted_p_value, upbm, tf, cosmic_tier,
           everything())
  
  data_reg_result_per_gene_top_expanded <- bind_rows(
    data_reg_result_per_gene_top_expanded,
    temp_result)
}

setwd(output.path)
filename_output <- paste0(datasetname, "_reg_result_weighted_top_gene_expanded.csv")
write.csv(data_reg_result_per_gene_top_expanded,
          filename_output, row.names = F, quote = F)

# output the gene each mutation is associated with
data_mutations_reg_annotated <- data_mutations_reg %>% 
  genome_inner_join(data_reg_annotated) %>%
  select(chromosome = chromosome.x,
         start = start.x,
         end = end.x,
         icgc_mutation_id, icgc_donor_id, ref, mut, gene) %>%
  distinct()

setwd(output.path)
filename_output <- paste0(datasetname, "_reg_mutations_annotated.csv")
write.csv(data_mutations_reg_annotated,
          filename_output, row.names = F, quote = F)
