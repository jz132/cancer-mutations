## Prepare genomic processing files

library(tidyverse)
library(data.table)
library(fuzzyjoin)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

rename <- dplyr::rename

# ============ PATH INFORMATION ============ 
source("functions.R")

promoter_path <- "dataset/epdata/promoters_refseq_coding.txt"
enhancer_path <- "dataset/epdata/all_enhancers_fantom.txt"
icgc_snps_path <- "dataset/icgc_processed/icgc_snps_LIRI_JP_no_exons.tsv"
datasetname <- "LIRI-JP"

out_tbl_dir <- "dataset/bg_entropy_data"
out_fig_dir <- "output_figures"

# change this depending on the input
cistype <- "enhancer" # enhancer or promoter

# ============ READ INPUT FILES ============ 
icgc_snps <- fread(icgc_snps_path)
num_donors <- nrow(icgc_snps %>% distinct(icgc_donor_id))

# fread automaticaly does read delim, see https://stackoverflow.com/questions/34144314/how-to-have-fread-perform-like-read-delim
promoters <- fread(promoter_path, col.names = c("chromosome", "start", "end", "promoter","strand")) %>%
  select(chromosome, start, end, promoter) %>% distinct()
enhancers <- fread(enhancer_path, col.names = c("chromosome", "start", "end", "enhancer"))

# ============ PROCESSING ============ 

cis <- if (cistype == "enhancer") enhancers else promoters
paste("number of all",cistype,":",nrow(cis))

# get genomic regions with mutations
cis_snps <- cis %>% 
  genome_left_join(icgc_snps, by = c("chromosome", "start", "end"))

# ============ PART1: GET MUTATIONS STATISTICS ============ 

# count how many mutations in the cis-regulatory regions
# length is inclusive
cis_snpcount <- cis_snps %>% 
    select(chromosome = chromosome.x,
                     start = start.x,
                     end = end.x,
                     icgc_mutation_id) %>%
    mutate(count = ifelse(is.na(icgc_mutation_id), 0, 1)) %>%
    group_by(chromosome, start, end) %>%
    summarise(count = sum(count)) %>%
    mutate(length = end - start + 1) %>%
    ungroup()
write.csv(cis_snpcount, paste0(out_tbl_dir,"/",cistype,"_snp_count.csv") , row.names = FALSE, quote=FALSE)
paste("Mutation count in",cistype)
num_cis_snpcount <- as.data.frame(cis_snpcount %>% group_by(count) %>% tally())
print(num_cis_snpcount, row.names = FALSE)

# Calculate the mutation rate 
snp_rate <- sum(cis_snpcount[["count"]]) / sum(cis_snpcount[["length"]]) / num_donors
paste("Mutation rate in",cistype,formatC(snp_rate, format = "e", digits = 2))

# Regulatory regions with SNPs
cis_with_snps <- cis_snpcount %>% filter(count > 0)
paste(cistype,"with mutations:",nrow(cis_with_snps))

# Get all the mutations that happen within the cis-regions
snps_in_cis <- cis_snps %>%
  na.omit() %>%
  select(chromosome = chromosome.y,
         start = start.y,
         end = end.y,
         icgc_mutation_id, icgc_donor_id,
         ref,mut) %>%
  distinct()
write.csv(snps_in_cis, paste0(out_tbl_dir,"/snps_in_",cistype,".csv") , row.names = FALSE, quote=FALSE)

# ============ PART2: GET TRINUCLEOTIDE FREQUENCY ============ 

genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg19")

# 1. FROM ALL CIS REGIONS, +/- 1 for the trinucleotide context
seqs_in_cis <- getSeq(genome, names = cis$chromosome,
                            start = cis$start - 1,
                            end = cis$end + 1)
cis_trimat <- reverseMerge(trinucleotideFrequency(seqs_in_cis))
cis_trifreq <- enframe(colSums(cis_trimat), name = "trinucleotide", value = "bg_count")

# 2. AROUND SNPS
# trinucleotide background of the SNPs
seqs_around_snps_ref <- getSeq(genome, names = snps_in_cis$chromosome,
                                     start = snps_in_cis$start - 1,
                                     end = snps_in_cis$end + 1)
seqs_around_snps_mut <- replaceLetterAt(seqs_around_snps_ref, 
                                         at = matrix(c(F, T, F), 
                                                     nrow = length(seqs_around_snps_ref),
                                                     ncol = 3,
                                                     byrow = T), 
                                         snps_in_cis$mut)
snps_trimat <- reverseMerge(trinucleotideFrequency(seqs_around_snps_ref))
snps_trifreq <- enframe(colSums(snps_trimat), name = "trinucleotide", value = "snp_count") 

# 3. COMBINING BACKGROUND AND SNPS TRINUCLEOTIDE

# ====== from ref to mut, and count ====== 
ref_mut_count <- tibble(
  ref = as.character(seqs_around_snps_ref),
  ref_rev = as.character(reverseComplement(seqs_around_snps_ref)),
  mut = as.character(seqs_around_snps_mut),
  mut_rev =  as.character(reverseComplement(seqs_around_snps_mut)),
) %>%
  mutate(ref = ifelse(ref < ref_rev, ref, ref_rev),
         mut = ifelse(ref < ref_rev, mut, mut_rev)) %>%
  group_by(ref, mut) %>%
  tally(name = "count") %>% 
  ungroup()

ggplot(data = ref_mut_count, mapping = aes(x = ref, y = count, 
                                           fill = paste0(substr(ref, 2, 2), "->", substr(mut, 2, 2)))) + 
  geom_col() + 
  labs(fill = "Mutation") +
  ggtitle(paste(cistype, " mutation")) + 
  theme(axis.text.x = element_text(angle = 45))
ggsave(paste0(out_fig_dir, "/", cistype,"_tri_mutdist.png"))
dev.off()

#  ====== Trinucleotide frequency and mutation rate  ====== 

# part where we calculate mut rate
trifreq_mut_rate <- cis_trifreq %>%
  inner_join(snps_trifreq, by="trinucleotide") %>%
  mutate(mut_rate = snp_count/(bg_count*num_donors))
write.csv(trifreq_mut_rate, paste0(out_tbl_dir,"/",cistype,"_tri_mutrate.csv") , row.names = FALSE, quote=FALSE)

ggplot(data = trifreq_mut_rate, aes(x = trinucleotide, y = bg_count)) + 
  geom_col() + 
  ggtitle(paste(cistype, "trinucleotide background distribution")) + 
  ylab("count") +
  theme(axis.text.x = element_text(angle = 45, size = 9), 
        axis.text.y = element_text(size = 10)
        )
ggsave(paste0(out_fig_dir,"/",cistype,"_tri_bg.png"))

dev.off()

ggplot(data = trifreq_mut_rate, aes(x = trinucleotide, y = mut_rate)) + 
  geom_col() + 
  ggtitle(paste("Trinucleotide mutation rate in ", cistype)) + 
  theme(axis.text.x = element_text(angle = 45, size = 9), 
        axis.text.y = element_text(size = 10)
  ) +
  labs(x = "", y = "mutation rate")
ggsave(paste0(out_fig_dir,"/",cistype,"_trimutrate.png"))
dev.off()

# ======== GET THE PROPORTION OF EACH TRINUCLEOTIDE ======== 

ref_mut_rate <- ref_mut_count %>% 
  group_by(ref) %>% 
  mutate(proportion = count/sum(count)) %>% 
  ungroup() %>% 
  inner_join(trifreq_mut_rate %>% select(ref = trinucleotide, mut_rate), by = "ref") %>%
  mutate(tri_mut_rate = proportion*mut_rate) %>%
  rename(ref_tri = ref, mut_tri = mut)
# mutate(mut_type = row_number())
print(head(data.table(ref_mut_rate)), row.names = FALSE)
write.csv(ref_mut_rate, paste0(out_tbl_dir,"/",cistype,"_mut_rate.csv") , row.names = FALSE, quote=FALSE)
