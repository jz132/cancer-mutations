### In this code, we use primarily the enhancer-tss association file
### from the Fantom project to build a pelinks_sameTSS file. 
### This new file gives for each TSS, all the promoter and enhancer regions 
### associatited with it.

# source packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("genomeIntervals")
library(genomeIntervals)
library(dplyr)
library(fuzzyjoin)
library(gtools)

# input path
home.path <- getwd()
fantom.path <- paste0(home.path, "/FANTOM")
refseq.path <- paste0(home.path, "/RefSeq")
processed.path <- paste0(home.path, "/Processed")

# output path
output.path <- home.path

# global parameters
options(stringsAsFactors = FALSE)
half_p <- 2000 # half promoter region length
half_tss <- 500 # enhancers within this distance to tss are removed
half_exon <- 200 # enhancers within this distance to exons are removed

## Part 1. Explore enhancer-promoter association and check the numbers in the paper
setwd(fantom.path)
data_enhancers_rob <- as_tibble(read.table("robust_enhancers.bed", 
                                           header = F, sep = '\t'))
data_enhancers_per <- as_tibble(read.table("permissive_enhancers.bed", 
                                           header = F, sep = '\t', skip = 1))
data_association_refseq <- as_tibble(read.table("enhancer_tss_associations_complete.bed", 
                                      header = F, sep = '\t'))

e_all_rob <- unique(data_enhancers_rob$V4) #all robust enhancers
e_all_per <- unique(data_enhancers_per$V4) #all permissive enhancers
all(e_all_rob %in% e_all_per) #robust enhancers are included in permissive ones

# check if enhancers overlap
overlap <- FALSE
for(chr in unique(data_enhancers_per$V1)){
  data_check <- data_enhancers_per[data_enhancers_per$V1 == chr, 1:3]
  data_check <- data_check[order(data_check$V2), ] #sorted data
  for(i in 1:(nrow(data_check)-1)){
    if(data_check$V3[i] > data_check$V2[i+1]){
      print("Overlapped")
      overlap <- TRUE
    }
  }
}
print(ifelse(overlap, "there is overlap between enhancers", "no overlap between enhancers"))

refseq_info_split <- strsplit(as.character(data_association_refseq$V4), split = ";")
e_eplinks_refseq <- unlist(lapply(refseq_info_split, `[[`, 1))
tss_eplinks_refseq <- unlist(lapply(refseq_info_split, `[[`, 2))
g_eplinks_refseq <- unlist(lapply(refseq_info_split, `[[`, 3))

all(e_eplinks_refseq %in% e_all_rob) #FALSE
all(e_eplinks_refseq %in% e_all_per) #TRUE 
# so eplinks_refseq is the enhancer-promoter links between 
# all permissive enhancers (not robust enhancers) and refseq TSS

# reproduce the numbers in the paper
length(tss_eplinks_refseq)/length(unique(e_eplinks_refseq)) #TSSs per enhancer
length(e_eplinks_refseq)/length(unique(tss_eplinks_refseq)) #enhancers per TSS
length(unique(e_eplinks_refseq))/length(e_all_per) #percentage of enhancers that have at least one TSS associated

## This is the end of using the original data ##
## It has several problems and is out of date ##

## Solved ##
# Because the original file contains missing data we need to clean the data before proceeding
# Already contacted FANTOM project staff and got this issue corrected
## Solved ##

## Part 2. Update the outdated data set and build promoter-enhancer links
setwd(refseq.path)
exons_region <- as_tibble(read.table("refseq_exons_171007.bed", sep = "\t")) %>%
  select(V1, V2, V3) %>%
  mutate(V2 = V2 - half_exon, V3 = V3 + half_exon) %>%
  rename(chromosome = V1, start = V2, end = V3)
tss_region <- as_tibble(read.table("refseq_TSS_hg19_170929.bed", sep = '\t')) %>%
  select(V1, V2, V3) %>%
  mutate(V2 = V2 - 500, V3 = V3 + 500) %>%
  rename(chromosome = V1, start = V2, end = V3)

data_enhancers_temp <- data_enhancers_per %>%
  select(V1, V2, V3) %>%
  rename(chromosome = V1, start = V2, end = V3)

data_enhancers_exclude <- data_enhancers_temp %>% 
  genome_inner_join(tss_region) %>%
  select(chromosome.x, start.x, end.x) %>%
  bind_rows(data_enhancers_temp %>% 
              genome_inner_join(exons_region) %>%
              select(chromosome.x, start.x, end.x)) %>%
  distinct() %>%
  rename(V1 = chromosome.x, V2 = start.x, V3 = end.x)

data_enhancers_new <- data_enhancers_per %>%
  anti_join(data_enhancers_exclude)
e_all_new <- unique(data_enhancers_new$V4)
all(e_eplinks_refseq %in% e_all_new) #need to exclude some enhancers from eplinks
list_drop <- which(!e_eplinks_refseq %in% e_all_new) #691 rows to exclude

# Update the enhancer, tss, and gene lists
refseq_info_split_filtered <- refseq_info_split[-list_drop]
e_eplinks_refseq <- unlist(lapply(refseq_info_split_filtered, `[[`, 1)) 
tss_eplinks_refseq <- unlist(lapply(refseq_info_split_filtered, `[[`, 2))
g_eplinks_refseq <- unlist(lapply(refseq_info_split_filtered, `[[`, 3))
tss_all <- unique(unlist(strsplit(tss_eplinks_refseq, split = ",")))

# Transform data from enhancer - promoters to promoter - enhancers
data_temp <- tibble(enhancer = e_eplinks_refseq, 
                    tss = tss_eplinks_refseq, 
                    gene = g_eplinks_refseq
)

mapping_tss_gene <- data_temp %>% 
  select(tss, gene) %>%
  distinct()

data_new <- data_temp %>% 
  select(tss, enhancer) %>%
  group_by(tss) %>%
  summarise(enhancer = paste(enhancer, collapse=";")) %>%
  inner_join(mapping_tss_gene, by = "tss") %>%
  distinct()

pelinks_raw <- data_new # 13588 promoters
colnames(pelinks_raw) <- c("tss", "enhancer", "gene_fantom")
pelinks_raw <- pelinks_raw %>% 
  mutate(e_count = sapply(strsplit(pelinks_raw$enhancer, ";"), length))
# these genes are annotated by fantom file enhancer_tss_associations.bed

# A histogram showing the number of enhancers each tss is associated with
setwd(output.path)
png(file = "ep_hist.png", width = 1200, height = 800, res = 160)
y <- hist(pelinks_raw$e_count)
plot(y, ylim=c(0, max(y$counts)+500), main = "Number of enhancers per RefSeq TSS",
     xlab = "Number of enhancers")
text(y$mids, y$counts+500, y$counts, cex=0.75)
dev.off()

# Read in the RefSeq data and Locate the TSSs
setwd(refseq.path)
data_tss <- as_tibble(read.table("refseq_TSS_hg19_170929.bed", 
                                 header = F, sep = '\t')) #all RefSeq tss
all(tss_all %in% data_tss$V4) #FALSE, 260 TSSs that no longer exist in the current RefSeq database

pelinks_sameTSS <- NULL
row_drop <- NULL
# check if all the transcripts have one and only one associated TSS
for(i in 1:nrow(pelinks_raw)){
  transcripts <- unlist(strsplit(pelinks_raw$tss[i], split = ","))
  base_start <- data_tss$V2[which(data_tss$V4 %in% transcripts)]
  strand <- as.character(data_tss$V6[which(data_tss$V4 %in% transcripts)])
  
  if(length(base_start) == 0){
    print(paste("RefSeq ID no longer valid", i))
    row_drop <- c(row_drop, i)
  } else if(!all(base_start == base_start[1])){
    print(paste("More than one TSS", i))
    row_drop <- c(row_drop, i)
  } else {
    pelinks_sameTSS <- rbind(pelinks_sameTSS, tibble(
      tss = pelinks_raw$tss[i],
      enhancer = pelinks_raw$enhancer[i],
      gene_fantom = pelinks_raw$gene_fantom[i],
      e_count = pelinks_raw$e_count[i],
      base_start = base_start[1],
      strand = strand[1]))
  }
}

length(row_drop) #363 rows with expired RefSeq ID or have more than one tss
pelinks_sameTSS <- pelinks_sameTSS %>%
  mutate(chr = gsub(":.*$", "", enhancer))

## Output raw promoters in bed format
p_raw <- data.frame(chr = pelinks_sameTSS$chr,
                    p_l = pelinks_sameTSS$base_start - half_p - 1,
                    p_u = pelinks_sameTSS$base_start + half_p,
                    tss = pelinks_sameTSS$tss)
p_raw <- as_tibble(p_raw)
filename_output <- paste0("promoters_raw_", half_p, ".bed")
setwd(output.path)
write.table(p_raw, filename_output, quote = F, sep = "\t",
            row.names = F, col.names = F)

## Use Galaxy to exclude the coding exons from each promoter and 
## get fragments of promoters

## Read in the promoters with coding exons excluded
setwd(processed.path)
p_noexons <- read.table("promoters_nocodingexons_2000_171102.bed", header = F, 
                        sep ='\t')
p_noexons <- as_tibble(p_noexons)
colnames(p_noexons) <- c("chr", "p_l", "p_u", "tss")

# Add the promoter variable in pelinks_sameTSS
pelinks_sameTSS <- pelinks_sameTSS %>% 
  inner_join(p_noexons) %>%
  mutate(promoter = paste0(chr, ":", p_l, "-", p_u)) %>%
  group_by(tss, enhancer, gene_fantom, e_count, base_start, strand, chr) %>%
  summarise(promoter = paste0(promoter, collapse = ";"))

promoters_associated <- subset(pelinks_sameTSS, select = c("tss", "promoter"))
enhancers_associated <- subset(pelinks_sameTSS, select = c("tss", "enhancer"))

# Output the promoters and enhancers that are assoicated with refseq TSSs
setwd(output.path)
p_temp <- unlist(strsplit(promoters_associated$promoter, ";"))
promoters_associated_u <- unique(p_temp)
# write.table(promoters_associated, "p_eplinks_refseq.txt",
#             row.names = F, col.names = F, quote = F)
write.table(promoters_associated_u, "p_eplinks_refseq_unique.txt",
            row.names = F, col.names = F, quote = F)

e_temp <- unlist(strsplit(enhancers_associated$enhancer, ";"))
enhancers_associated_u <- unique(e_temp)
# write.table(enhancers_associated, "e_eplinks_refseq.txt",
#             row.names = F, col.names = F, quote = F)
write.table(enhancers_associated_u, "e_eplinks_refseq_unique.txt",
            row.names = F, col.names = F, quote = F)

## Read in coding exons
setwd(refseq.path)
refseq_codingexons <- read.table("refseq_codingexons_171007.bed", sep = "\t")
refseq_codingexons$V4 <- gsub("_cds.*$", "", refseq_codingexons$V4)
for(i in 1:nrow(pelinks_sameTSS)){
  transcripts <- unlist(strsplit(pelinks_sameTSS$tss[i], split = ";"))
  refseq_codingexons_subset <- 
    refseq_codingexons[which(refseq_codingexons$V4 %in% transcripts),]
  cds <- paste0(refseq_codingexons_subset$V1, ":",
                refseq_codingexons_subset$V2, "-",
                refseq_codingexons_subset$V3,
                collapse = ";")
  pelinks_sameTSS$cds[i] <- cds
}

pelinks_sameTSS <- pelinks_sameTSS %>% 
  select("tss", "chr", "base_start", "promoter", "enhancer", 
         "cds", "gene_fantom", "strand", "e_count")

## Output the result: TSS, promoter, enhancer, coding exons, etc.
setwd(output.path)
write.table(pelinks_sameTSS, "pelinks_sameTSS.txt", quote = F,
            row.names = F)
write.csv(pelinks_sameTSS, "pelinks_sameTSS.csv", quote = F,
          row.names = F)
# Note: this file contains NA in some promoters


# Try to use the promoters annotated by FANTOM itself
# data_temp <- data_association_cell
# p_eplinks_cell <- data_association_cell$promoter
# p_eplinks_cell_u <- unique(p_eplinks_cell)
# data_new <- NULL
# for(p in p_eplinks_cell_u[1:1000]){
#   row_select <- which(data_temp$promoter == p)
#   enhancer_associated <- paste(data_temp$enhancer[row_select], collapse = ";")
#   data_new <- rbind(data_new, c(as.character(p), enhancer_associated))
#   data_temp <- data_temp[-row_select, ]
# }
# colnames(data_new) <- c("promoter", "enhancer")
# data_new <- data.frame(data_new)
# data_new$promoter <- as.character((data_new$promoter))
# data_new$enhancer <- as.character((data_new$enhancer))


