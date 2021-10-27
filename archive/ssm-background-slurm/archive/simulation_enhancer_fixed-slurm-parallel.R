source("~/r_projects/cancer-mutations/ssm-background-slurm/trinucleotide_bg_enhancer-slurm.R")
print("successfully sourced trinucleotide_bg_enhancer-slurm.R")

slice <- dplyr::slice

## arguments from the bash file ##
args <- commandArgs(trailingOnly = TRUE)
cmd_arg <- as.numeric(args[1])

## functions ##
## absmax() gets the element with maximum absolute value from a vector of numbers
absmax <- function(x) {x[which.max(abs(x))]}

## DNAToBin() transforms DNA to binary number
## can apply R function strtoi() to further transform binary to decimal number
DNAToBin <- function(DNA){
  temp <- DNA
  temp <- gsub("A", "00", temp)
  temp <- gsub("C", "01", temp)
  temp <- gsub("G", "10", temp)
  temp <- gsub("T", "11", temp)
  return(temp)
}

# file paths and file names
synthetic.mutation.path <- "/data/gordanlab/jingkang/Desktop/Gordanlab/Data/synthetic_mutations"
qbic.model.path <- "/data/gordanlab/jingkang/Desktop/Gordanlab/Data/qbic/predmodel"
output.path <- "~/r_projects/cancer-mutations/ssm-background-slurm/output"

## Part 1: Import all possible single mutations for enhancers with mutations ##
enhancers_w_mutation <- which(data_enhancers_mutated$count > 0)

setwd(synthetic.mutation.path)
data_enhancer_mutation_sim <- read_tsv("enhancer_mutation_sim.vcf", 
                                       col_names = c("chromosome", "pos", "enhancer", "ref", "mut"))
seq_enhancer_mutation_sim <- read_table("seq_enhancer_mutation_sim.txt", col_names = F) %>% pull()

## Part 2: Predict the effect of the generated single mutations using QBiC ##
setwd(qbic.model.path)
table_12mer_myc <- read.table("prediction6mer.Homo_sapiens|NA|Shen2018|MYC-MAX.txt", 
                              header = T, fill = T, blank.lines.skip = F)
table_12mer_myc <- as_tibble(table_12mer_myc) %>% 
  mutate(idx = row_number() - 1)

print("table_12mer_myc imported")

data_enhancer_mutation_sim <- data_enhancer_mutation_sim %>% 
  mutate(twimer = paste0(seq_enhancer_mutation_sim, mut)) %>%
  mutate(idx = DNAToBin(twimer)) %>%
  mutate(idx = strtoi(idx, base = 2))

print("data_enhancer_mutation_sim processed")

data_enhancer_mutation_prediction <- data_enhancer_mutation_sim %>% 
  inner_join(table_12mer_myc) %>%
  mutate(ref_tri = substr(twimer, 5, 7)) %>%
  mutate(mut_tri = paste0(substr(ref_tri, 1, 1), mut, substr(ref_tri, 3, 3))) %>% 
  mutate(ref_tri_rev = as.character(reverseComplement(DNAStringSet(ref_tri))),
         mut_tri_rev = as.character(reverseComplement(DNAStringSet(mut_tri)))) %>% 
  mutate(ref_tri = ifelse(ref_tri < ref_tri_rev, ref_tri, ref_tri_rev),
         mut_tri = ifelse(ref_tri < ref_tri_rev, mut_tri, mut_tri_rev)) %>% 
  select(-c(ref_tri_rev, mut_tri_rev))

print("key join completed")

## Part 3: Simulation ##
set.seed(1234*cmd_arg)
n_run <- 10^4

# define mutation type based on the alphabetical order of the trinucleotides
table_mutation_tri_mut_type <- table_mutation_tri %>%
  select(ref_tri = ref, 
         mut_tri = mut) %>% 
  mutate(mut_type = row_number()) 

# add this mutation type information in our prediction table
data_enhancer_mutation_prediction <- data_enhancer_mutation_prediction %>% 
  inner_join(table_mutation_tri_mut_type)

# mutation rate for each of the 32 trinucleotides
tri_mut_rate <- freq_tri %>% pull(mut_rate) 
names(tri_mut_rate) <- freq_tri %>% pull(trinucleotide)

# proportion of 3 possible mutations for each of the 32 trinucleotides
mat_prop <- matrix(table_mutation_tri %>% 
                     group_by(ref) %>% 
                     mutate(prop = count/sum(count)) %>% 
                     pull(prop),
                   nrow = 32, ncol = 3, byrow = T) 
rownames(mat_prop) <- names(tri_mut_rate)

p_value_list <- rep(1, nrow(data_enhancers_mutated))
# which.max(data_enhancers_mutated$count)

for(enhancer_ex in enhancers_w_mutation){
  print(enhancer_ex)
  
  # the count of each trinucleotide in the enhancer
  vec_tri <- reverseMerge(trinucleotideFrequency(seq_enhancers_fantom[enhancer_ex]))
  n_tri <- length(vec_tri)
  
  # first random sampling from a multinomial distribution for 
  # the number of mutated trinucleotide, controlling the same number of mutations as real patients
  n_mut_sim <- t(rmultinom(n_run, data_enhancers_mutated$count[enhancer_ex], vec_tri*tri_mut_rate))
  colnames(n_mut_sim) <- names(vec_tri)
  
  # then random sampling from a multinomial distribution for
  # each trinucleotide to trinucleotide mutation
  n_mut_sim_expanded <- matrix(0, nrow = n_run, ncol = 3*n_tri)
  for(i in 1:n_run){
    for(j in 1:ncol(n_mut_sim)){
      if(n_mut_sim[i,j] != 0){
        n_mut_sim_expanded[i,(3*j-2):(3*j)] <- rmultinom(1, n_mut_sim[i,j], mat_prop[j,])
      }
    }
  }
  
  # a summary table for the number of trinucleotide to trinucleotide mutations
  table_mutation_tri_result <- table_mutation_tri %>%
    select(ref, mut) %>%
    mutate(sim_count = colSums(n_mut_sim_expanded))
  
  # mutation effect prediction for the simulated results
  data_prediction_ex <- data_enhancer_mutation_prediction %>%
    filter(enhancer == enhancer_ex)
  mut_sim_effect_table <- n_mut_sim_expanded
  for(i in 1:nrow(n_mut_sim_expanded)){
    for(j in 1:ncol(n_mut_sim_expanded)){
      n_mut <- n_mut_sim_expanded[i,j]
      if(n_mut == 0){
        mut_sim_effect_table[i,j] <- 0
      }
      else{
        data_temp <- data_prediction_ex %>%
          filter(mut_type == j)
        if(nrow(data_temp) != 0){
          sim_result <- sample(data_temp$diff, n_mut, replace = T)
          mut_sim_effect_table[i,j] <- absmax(sim_result)
        }
      }
    }
  }
  
  data_enhancers_mutated_ex <- data_enhancers_mutated[enhancer_ex, ]
  
  mut_enhancer_actual <- mut_enhancer %>%
    genome_inner_join(data_enhancers_mutated_ex) %>%
    select(chromosome = chromosome.x,
           pos = start.x,
           icgc_mutation_id,
           icgc_donor_id,
           ref,
           mut) %>%
    inner_join(data_prediction_ex)
  
  mut_enhancer_actual_effect <- absmax(mut_enhancer_actual$diff)
  mut_enhancer_sim_effect <- apply(mut_sim_effect_table, 1, absmax)
  p_value <- sum(abs(mut_enhancer_sim_effect) >= abs(mut_enhancer_actual_effect))/n_run
  p_value_list[enhancer_ex] <- p_value
  
  print(p_value)
}

setwd(output.path)
write.table(p_value_list,  paste0("p_value_list_enhancer_fixed_", cmd_arg, ".txt"), 
            row.names = F, col.names = F)

