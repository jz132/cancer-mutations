## Subset the input dataset with the information we need

library(dplyr)
library(data.table)

# ============ PATH INFORMATION ============ 
setwd("/Users/vincentiusmartin/Research/CancerMutation")
icgc_mut_path <- "/Volumes/VINCENTIUS_TOSHIBA/genomic_data/CancerMutation/simple_somatic_mutation.open.LIRI-JP.tsv"
icgc_mut_outpath <- "dataset/input/icgc_snps_LIRI_JP"

icgc_exp_path <- "/Volumes/VINCENTIUS_TOSHIBA/genomic_data/CancerMutation/exp_seq.LIRI-JP.tsv"
icgc_exp_outpath <- "dataset/input/icgc_exp_LIRI_JP"

exon_path <- "dataset/input/all_exons_refseq.txt"

# set chromosome order for sorting
chrom_order<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")

# ============ READ THE SNPS DATA ============ 
# we only focus on SNPs
icgc_mut <- fread(icgc_mut_path, select=c("chromosome","chromosome_start","chromosome_end",
                                  "icgc_mutation_id", "icgc_donor_id","reference_genome_allele",
                                  "mutated_from_allele","mutated_to_allele","mutation_type",
                                  "sequencing_strategy","consequence_type")) %>%
        filter(sequencing_strategy == "WGS", mutation_type == "single base substitution", !is.na(consequence_type)) %>%
        select(chromosome,
               start = chromosome_start,
               end = chromosome_end,
               icgc_mutation_id,
               icgc_donor_id,
               ref = reference_genome_allele,
               mutated_from_allele,
               mut = mutated_to_allele) %>%
        distinct() %>%
        mutate(chromosome = paste("chr",chromosome,sep=""))
if(!all(icgc_mut$mutated_from_allele == icgc_mut$ref)){
  warning("reference genome allele is not the same as mutated from allele")
}
icgc_mut <- icgc_mut %>% select(-mutated_from_allele)
icgc_mut$chromosome <- factor(icgc_mut$chromosome, levels=chrom_order)
icgc_mut <- icgc_mut %>% arrange(icgc_donor_id,chromosome,start) 

cat(sprintf("Nrow icgc mut: %d\n", nrow(icgc_mut)))
cat(sprintf("Num donors: %d\n", length(unique(icgc_mut$icgc_donor_id))))
cat(sprintf("Num mutations: %d\n", length(unique(icgc_mut$icgc_mutation_id))))
write.table(icgc_mut, paste0(icgc_mut_outpath,".tsv"), quote=FALSE, row.names=FALSE,sep="\t")

# filter out mutations within exons
exons <- fread(exon_path, col.names = c("chromosome", "start", "end", "exons")) %>% select(chromosome, start, end)
icgc_mut_no_exons <- icgc_mut %>%
  genome_anti_join(exons, by = c("chromosome", "start", "end"))
write.table(icgc_mut_no_exons, paste0(icgc_mut_outpath,"_no_exons.tsv"), quote=FALSE, row.names=FALSE,sep="\t")

# ============ READ THE RNA-SEQ DATA ============ 

wgs_donors <- unique(icgc_mut$icgc_donor_id)
icgc_exp <- fread(icgc_exp_path, 
                  select=c("icgc_donor_id", "gene_id","normalized_read_count", "analysis_id")) %>%
            distinct()
cat(sprintf("Num donors from icgc exp: %d\n", length(unique(icgc_exp$icgc_donor_id))))
            
icgc_exp <- icgc_exp %>% filter(icgc_donor_id %in% wgs_donors)
cat(sprintf("Overlapping icgc donor from icgc exp: %d\n", length(unique(icgc_exp$icgc_donor_id))))
cat(sprintf("Nrow icgc exp with donor overlap: %d\n", nrow(icgc_exp)))
write.table(icgc_exp, icgc_exp_outpath, quote=FALSE, row.names=FALSE, sep="\t")


