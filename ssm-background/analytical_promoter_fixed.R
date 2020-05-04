t0 <- proc.time()

source("~/r_projects/cancer-mutations/ssm-background/trinucleotide_bg_promoter.R")
rm(list = ls(pattern = "^data_icgc")) # free some space

t1 <- proc.time()
print(t1 - t0)

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
filename_table_12mer <- "prediction6mer.Homo_sapiens|NA|Shen2018|MYC-MAX.txt"

# import the 12-mer prediction table
setwd(qbic.model.path)
# wondering if there is a way to use read_table() to correctly read this table in
table_12mer <- read.table(filename_table_12mer, header = T, fill = T, 
                          blank.lines.skip = F)
table_12mer <- as_tibble(table_12mer) %>% 
  mutate(idx = row_number() - 1)

# proportion of 3 possible mutations for each of the 32 trinucleotides
table_mutation_tri_mut_rate <- table_mutation_tri %>% 
  group_by(ref) %>% 
  mutate(prop = count/sum(count)) %>% 
  ungroup() %>% 
  inner_join(freq_tri %>% select(ref = trinucleotide, mut_rate)) %>%
  mutate(tri_mut_rate = prop*mut_rate) %>% 
  select(ref_tri = ref,
         mut_tri = mut,
         tri_mut_rate) %>% 
  mutate(mut_type = row_number())

promoters_w_mutation <- which(data_promoters_mutated$count > 0)
mut_effect <- rep(0, nrow(data_promoters_mutated))
p_value_list <- rep(1, nrow(data_promoters_mutated))

t2 <- proc.time()
print(t2 - t1)

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
    mutate(diff = table_12mer$diff[idx+1]) %>%
    mutate(ref_tri = substr(twimer, 5, 7)) %>%
    mutate(mut_tri = paste0(substr(ref_tri, 1, 1), mut, substr(ref_tri, 3, 3))) %>% 
    mutate(ref_tri_rev = as.character(reverseComplement(DNAStringSet(ref_tri))),
           mut_tri_rev = as.character(reverseComplement(DNAStringSet(mut_tri)))) %>% 
    mutate(ref_tri = ifelse(ref_tri < ref_tri_rev, ref_tri, ref_tri_rev),
           mut_tri = ifelse(ref_tri < ref_tri_rev, mut_tri, mut_tri_rev)) %>% 
    select(-c(ref_tri_rev, mut_tri_rev))  %>% 
    inner_join(table_mutation_tri_mut_rate) %>% 
    mutate(cond_tri_mut_rate = tri_mut_rate/sum(tri_mut_rate))
  
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
  n_mutations_actual <- nrow(mut_promoter_actual)
  
  # a subset of all possible promoter mutations with less effect than actually observed
  data_promoter_less_than_actual <- data_promoter_mutation_prediction %>%
    filter(abs(diff) < abs(mut_promoter_actual_effect))
  p_value <- 1 - sum(data_promoter_less_than_actual$cond_tri_mut_rate)^n_mutations_actual
  
  mut_effect[promoter_ex] <- mut_promoter_actual_effect
  p_value_list[promoter_ex] <- p_value
  
  print(mut_promoter_actual_effect)
  print(p_value)
}

t3 <- proc.time()
print(t3 - t2)

setwd(output.path)
write.table(p_value_list, paste0("p_value_promoter_", filename_table_12mer),
            row.names = F, col.names = F)

write.table(mut_effect, paste0("mut_effect_promoter_", filename_table_12mer),
            row.names = F, col.names = F)
