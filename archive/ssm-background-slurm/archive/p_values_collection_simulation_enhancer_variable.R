library(gtools)
setwd("/Users/jz132/Desktop/hardac-xfer/p_value_list_variable/enhancer/")
filenames <- mixedsort(Sys.glob("*.txt"))

temp_filename <- filenames[1]
n_enhancers <- nrow(read.table(temp_filename))

p_value_list <- numeric(n_enhancers)
for(filename in filenames){
  new_p_value <- read.table(filename)[,1]
  p_value_list <- p_value_list + new_p_value/length(filenames)
}

p_values <- (p_value_list*10^6 + 1)/(10^6 + 1) # correct the monte carlo p-values
p_values_non_one <- p_values[p_values < 1]
hist(p_values)
hist(p_values_non_one)

x <- qunif(ppoints(length(p_values)))
qqplot(-log(x, base = 10), -log(p_values, base = 10), 
       main = "QQ-plot of Enhancer P-values", 
       xlab = "expected -log(p, base=10)",
       ylab = "observed -log(p, base=10)")
abline(a = 0, b = 1)

data_enhancers_mutated <- data_enhancers_mutated %>% 
  mutate(p_value = p_values) %>%
  mutate(log.p_value = -log(p_value, 10)) %>% 
  arrange(desc(log.p_value)) %>%
  mutate(expected_logp = -log(x, base = 10))

setwd("/Users/jz132/r_projects/cancer-mutations/pelinks")
data_eplinks_short <- read_delim("eplinks-fantom-filtered.csv", delim = ",")
data_eplinks_long <- data_eplinks_short %>%
  mutate(tss = strsplit(tss, ";")) %>%
  unnest(tss)

