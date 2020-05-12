source("~/r_projects/cancer-mutations/ssm-background/trinucleotide_bg_promoter.R")
library(gtools)

tf_name <- "TCF7L2"
input.path <- paste0("/Users/jz132/Desktop/hardac-xfer/",
                     tf_name,
                     "/p_value_list_fixed_analytical/promoter/")

# collect the p-values result
setwd(input.path)
filenames <- mixedsort(Sys.glob("*.txt"))

n_promoters <- nrow(data_promoters_refseq)
promoters_w_mutation <- which(data_promoters_mutated$count > 0)

p_value_list <- rep(1, n_promoters)
for(filename in filenames){
  new_p_value <- read_delim(filename, delim = " ", col_names = "p_value") %>% 
    pull(p_value)
  p_value_list <- pmin(p_value_list, new_p_value)
}

p_values_analytical <- p_value_list[promoters_w_mutation]

x <- qunif(ppoints(length(p_values_analytical)))
qqplot(-log(x, base = 10), -log(p_values_analytical, base = 10), 
       main = "QQ-plot of promoter P-values", 
       xlab = "expected -log(p, base=10)",
       ylab = "observed -log(p, base=10)")
abline(a = 0, b = 1)

data_promoters_mutated_effect <- data_promoters_mutated %>% 
  mutate(p_value = p_value_list) %>% 
  arrange(p_value)

top_promoters <- data_promoters_mutated_effect %>% 
  select(chromosome, start, end, promoter) %>% 
  mutate(pos = end - 1000) %>%
  slice(1:20)

data_promoters_mutated_effect_output <- data_promoters_mutated_effect %>% 
  filter(count > 0)

write.csv(data_promoters_mutated_effect_output,
          file = paste0(tf_name, "_binding_promoter_p_value.csv"),
          row.names = F)
