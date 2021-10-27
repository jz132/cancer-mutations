source("~/r_projects/cancer-mutations/ssm-background-slurm/trinucleotide_bg_promoter-slurm.R")
print("successfully sourced trinucleotide_bg_promoter-slurm.R")

slice <- dplyr::slice
rename <- dplyr::rename

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
qbic.model.path <- "/data/gordanlab/vincentius/qbic/predmodel"
output.path <- "~/r_projects/cancer-mutations/ssm-background-slurm/output"

# import the 12-mer prediction table
setwd(qbic.model.path)
filename_table_12mer <- Sys.glob("*.txt")[cmd_arg]

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

result_list <- list()
for(i in 1:length(promoters_w_mutation)){
  promoter_ex <- promoters_w_mutation[i]
  print(promoter_ex)
  
  ## Part 1: Generate all possible single mutations for promoters with mutations ##
  # get the location and sequence of the promoter
  data_promoter_ex <- data_promoters_refseq %>% slice(promoter_ex) # the promoter location
  seq_promoter_ex <- seq_promoters_refseq[[promoter_ex]] # the promoter sequence
  seq_promoter_ex_bg <- getSeq(genome, data_promoter_ex$chromosome, 
                               data_promoter_ex$start - 5, 
                               data_promoter_ex$end + 5) # the promoter sequence including 5bp context
  seq_promoter_ex_11mer <- DNAStringSet(seq_promoter_ex_bg, 
                                        start = 1:(length(seq_promoter_ex_bg)-10),
                                        width = 11) # cut the promoter into 11mers
  
  # generate all possible single mutations
  data_promoter_all_possible_mutations <- tibble(
    chromosome = data_promoter_ex$chromosome,
    pos = seq(data_promoter_ex$start, data_promoter_ex$end),
    promoter = promoter_ex,
    seq_mut_bg = as.character(seq_promoter_ex_11mer)
  ) %>% 
    mutate(ref = unlist(strsplit(as.character(seq_promoter_ex), ""))) %>% 
    slice(rep(seq(1, length(seq_promoter_ex)), each = 4)) %>% 
    mutate(mut = rep(c("A", "C", "G", "T"), length(seq_promoter_ex))) %>% 
    filter(ref != mut) %>% 
    mutate(twimer = paste0(seq_mut_bg, mut)) %>%
    mutate(idx = DNAToBin(twimer)) %>%
    mutate(idx = strtoi(idx, base = 2)) %>%
    select(-seq_mut_bg)
  
  ## Part 2: Predict the effect of the synthetic mutations using QBiC ##
  # use the 12-mer table to predict the effect
  data_promoter_mutation_prediction <- data_promoter_all_possible_mutations %>% 
    mutate(effect = table_12mer$diff[idx+1]) %>%
    mutate(ref_tri = substr(twimer, 5, 7)) %>%
    mutate(mut_tri = paste0(substr(ref_tri, 1, 1), mut, substr(ref_tri, 3, 3))) %>% 
    mutate(ref_tri_rev = as.character(reverseComplement(DNAStringSet(ref_tri))),
           mut_tri_rev = as.character(reverseComplement(DNAStringSet(mut_tri)))) %>% 
    mutate(ref_tri = ifelse(ref_tri < ref_tri_rev, ref_tri, ref_tri_rev),
           mut_tri = ifelse(ref_tri < ref_tri_rev, mut_tri, mut_tri_rev)) %>% 
    select(-c(ref_tri_rev, mut_tri_rev))  %>% 
    inner_join(table_mutation_tri_mut_rate, by = c("ref_tri", "mut_tri")) %>% 
    mutate(cond_tri_mut_rate = tri_mut_rate/sum(tri_mut_rate))
  
  # mutation effect prediction for the actual mutations in the promoter
  mut_promoter_effect_prediction <- mut_promoter %>% 
    genome_inner_join(data_promoter_ex,
                      by = c("chromosome", "start", "end")) %>%
    select(chromosome = chromosome.x,
           pos = start.x,
           icgc_mutation_id,
           icgc_donor_id,
           ref,
           mut) %>%
    inner_join(data_promoter_mutation_prediction,
               by = c("chromosome", "pos", "ref", "mut")) %>%
    mutate(p_less = 0)
  
  for(j in 1:nrow(mut_promoter_effect_prediction)){
    mut_promoter_effect_prediction$p_less[j] <- sum(
      data_promoter_mutation_prediction %>%
        filter(effect < mut_promoter_effect_prediction$effect[j]) %>%
        pull(cond_tri_mut_rate)
    )
  }
  
  result <- mut_promoter_effect_prediction %>% 
    select(icgc_mutation_id, icgc_donor_id, promoter, effect, p_less)
  
  result_list[[i]] <- result
}

promoter_result <- bind_rows(result_list)

setwd(output.path)
write.table(promoter_result, paste0("promoter_mutation_result_", filename_table_12mer),
            row.names = F)
