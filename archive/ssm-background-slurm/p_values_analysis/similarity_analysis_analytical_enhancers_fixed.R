# analysis on whether different tables predict similarly for corresponding tf
source("~/r_projects/cancer-mutations/ssm-background/trinucleotide_bg_enhancer.R")
library(psych)

setwd("/Users/jz132/Desktop/hardac-xfer/mappings")
mapping_table_per_tf <- read_csv("mapping_to_12mer_tbl.csv") %>% 
  select(tf = gene, 
         upbm) %>% 
  mutate(upbm = gsub(" ", "", upbm, fixed = T)) %>% 
  mutate(upbm = gsub("prediction6mer.Mus_musculus|M01338_1.94d|Badis09|Six6_2267.4=v1.txt",
                     "prediction6mer.Mus_musculus|M01339_1.94d|Berger08|Six6_2267.4.txt",
                     upbm, fixed = T))

tf_names <- mapping_table_per_tf$tf
for(tf_name in tf_names){
  upbm <- mapping_table_per_tf %>%
    filter(tf == tf_name) %>%
    pull(upbm)
  
  upbm_filenames <- unlist(strsplit(upbm, split = ";"))
  p_value_filenames <- paste0("p_value_enhancer_", upbm_filenames)
  mut_effect_filenames <- paste0("mut_effect_enhancer_", upbm_filenames)
  num_files <- length(p_value_filenames)
  
  data_enhancers_mutated_effect_list <- list()
  setwd("/Users/jz132/Desktop/hardac-xfer/raw_data/enhancer")
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
    
    data_enhancers_mutated_effect_list[[i]] <- data_enhancers_mutated_effect
  }
  
  names(data_enhancers_mutated_effect_list) <- upbm_filenames
  print(data_enhancers_mutated_effect_list)
  
  cor_plot_names <- unlist(lapply(strsplit(upbm_filenames, split = "|", fixed = T), "[", 4))
  
  if(length(data_enhancers_mutated_effect_list) > 1){
    table_mut_effect <- matrix(0, nrow = length(enhancers_w_mutation),
                               ncol = length(data_enhancers_mutated_effect_list))
    for(j in 1:length(data_enhancers_mutated_effect_list)){
      table_mut_effect[,j] <- data_enhancers_mutated_effect_list[[j]]$mut_effect
    }
    colnames(table_mut_effect) <- cor_plot_names
    setwd("/Users/jz132/Desktop/hardac-xfer/plots")
    png(file = paste0(tf_name, "_enhancer.png"), width = 1600, height = 1600, res = 160)
    pairs.panels(table_mut_effect, digits = 4, smooth = F, ellipses = F,
                 main = tf_name)
    dev.off()
  }
}
