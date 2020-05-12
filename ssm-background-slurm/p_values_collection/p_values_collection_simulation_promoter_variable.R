source("~/r_projects/cancer-mutations/ssm-background/trinucleotide_bg_promoter.R")

library(gtools)
setwd("/Users/jz132/Desktop/hardac-xfer/p_value_list_variable_simulation/promoter/")
filenames <- mixedsort(Sys.glob("*.txt"))

n_promoters <- nrow(data_promoters_refseq)

p_value_list <- rep(1, n_promoters)
for(filename in filenames){
  new_p_value <- read_delim(filename, delim = " ", col_names = "p_value") %>% 
    pull(p_value)
  p_value_list <- pmin(p_value_list, new_p_value)
}

p_values <- (p_value_list*10^6 + 1)/(10^6 + 1)

x <- qunif(ppoints(length(p_values)))
qqplot(-log(x, base = 10), -log(p_values, base = 10), 
       main = "QQ-plot of Promoter P-values", 
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

