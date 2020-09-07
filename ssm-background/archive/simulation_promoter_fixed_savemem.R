source("~/r_projects/cancer-mutations/ssm-background/trinucleotide_bg_promoter.R")
rm(list = ls(pattern = "^data_icgc")) # free some space

options(tibble.width = Inf)
slice <- dplyr::slice

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
qbic.model.path <- "~/Desktop/Gordanlab/Data/qbic/predmodel"
output.path <- "~/r_projects/cancer-mutations/ssm-background/output"

# import the 12-mer prediction table
setwd(qbic.model.path)
# wondering if there is a way to use read_table() to correctly read this table in
table_12mer <- read.table("prediction6mer.Homo_sapiens|NA|Shen2018|MYC-MAX.txt", 
                              header = T, fill = T, blank.lines.skip = F)
table_12mer <- as_tibble(table_12mer) %>% 
  mutate(idx = row_number() - 1)

# define mutation type based on the alphabetical order of the trinucleotides
table_mutation_tri_mut_type <- table_mutation_tri %>%
  select(ref_tri = ref, 
         mut_tri = mut) %>% 
  mutate(mut_type = row_number())

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

## simulation code ##
set.seed(1234)
n_run <- 10^2
promoters_w_mutation <- which(data_promoters_mutated$count > 0)
p_value_list <- rep(1, nrow(data_promoters_mutated))

for(promoter_ex in promoters_w_mutation){
  print(promoter_ex)
  
  ## Part 1: Generate all possible single mutations for promoters with mutations ##
  # get the location and sequence of the promoter
  data_promoter_ex <- data_promoters_refseq %>% slice(promoter_ex) # the promoter location
  seq_promoter_ex <- seq_promoters_refseq[[promoter_ex]] # the promoter sequence
  
  # generate all possible mutations and store in a vcf format table
  data_promoter_mutation_ex_vcf <- tibble(
    chromosome = data_promoter_ex$chromosome,
    pos = seq(data_promoter_ex$start, data_promoter_ex$end),
    promoter = promoter_ex
  ) %>% 
    mutate(ref = unlist(strsplit(as.character(seq_promoter_ex), ""))) %>% 
    slice(rep(seq(1, length(seq_promoter_ex)), each = 4)) %>% 
    mutate(mut = rep(c("A", "C", "G", "T"), length(seq_promoter_ex))) %>% 
    filter(ref != mut)
  
  # get the 11-mer mutation context for all synthetic mutations
  seq_promoter_mutation_sim <- getSeq(genome,
                                      names = data_promoter_mutation_ex_vcf$chromosome, 
                                      start = data_promoter_mutation_ex_vcf$pos - 5,
                                      end = data_promoter_mutation_ex_vcf$pos + 5)
  
  ## Part 2: Predict the effect of the synthetic mutations using QBiC ##
  data_promoter_mutation_sim <- data_promoter_mutation_ex_vcf %>% 
    mutate(twimer = paste0(seq_promoter_mutation_sim, mut)) %>%
    mutate(idx = DNAToBin(twimer)) %>%
    mutate(idx = strtoi(idx, base = 2))
  
  # use the 12-mer table to predict the effect
  data_promoter_mutation_prediction <- data_promoter_mutation_sim %>% 
    inner_join(table_12mer) %>%
    mutate(ref_tri = substr(twimer, 5, 7)) %>%
    mutate(mut_tri = paste0(substr(ref_tri, 1, 1), mut, substr(ref_tri, 3, 3))) %>% 
    mutate(ref_tri_rev = as.character(reverseComplement(DNAStringSet(ref_tri))),
           mut_tri_rev = as.character(reverseComplement(DNAStringSet(mut_tri)))) %>% 
    mutate(ref_tri = ifelse(ref_tri < ref_tri_rev, ref_tri, ref_tri_rev),
           mut_tri = ifelse(ref_tri < ref_tri_rev, mut_tri, mut_tri_rev)) %>% 
    select(-c(ref_tri_rev, mut_tri_rev))  %>% 
    inner_join(table_mutation_tri_mut_type)
  
  ## Part 3: Simulation ##
  # the count of each trinucleotide in the promoter
  vec_tri <- reverseMerge(trinucleotideFrequency(seq_promoters_refseq[promoter_ex]))
  n_tri <- length(vec_tri)
  
  # first random sampling from a multinomial distribution for 
  # the number of mutated trinucleotide, controlling the same number of mutations as real patients
  n_mut_sim <- t(rmultinom(n_run, data_promoters_mutated$count[promoter_ex], vec_tri*tri_mut_rate))
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
  # table_mutation_tri_result <- table_mutation_tri %>% 
  #   select(ref, mut) %>%
  #   mutate(sim_count = colSums(n_mut_sim_expanded))
  
  # mutation effect prediction for the simulated results
  mut_sim_effect_table <- n_mut_sim_expanded
  for(i in 1:nrow(n_mut_sim_expanded)){
    for(j in 1:ncol(n_mut_sim_expanded)){
      n_mut <- n_mut_sim_expanded[i,j]
      if(n_mut == 0){
        mut_sim_effect_table[i,j] <- 0
      }
      else{
        data_temp <- data_promoter_mutation_prediction %>%
          filter(mut_type == j)
        if(nrow(data_temp) != 0){
          sim_result <- sample(data_temp$diff, n_mut, replace = T)
          mut_sim_effect_table[i,j] <- absmax(sim_result)
        }
      }
    }
  }
  
  # mutation effect prediction for the actual data
  mut_promoter_actual <- mut_promoter %>% 
    genome_inner_join(data_promoter_ex) %>%
    select(chromosome = chromosome.x,
           pos = start.x,
           icgc_mutation_id,
           icgc_donor_id,
           ref,
           mut) %>%
    inner_join(data_promoter_mutation_prediction)
  
  mut_promoter_actual_effect <- absmax(mut_promoter_actual$diff) 
  mut_promoter_sim_effect <- apply(mut_sim_effect_table, 1, absmax)
  p_value <- sum(abs(mut_promoter_sim_effect) >= abs(mut_promoter_actual_effect))/n_run
  p_value_list[promoter_ex] <- p_value
  
  print(p_value)
}

# setwd(output.path)
# write.table(p_value_list, "p_value_list_promoter_fixed_savemem.txt", row.names = F, col.names = F)
