library(gtools)

# collect the analytical p-values result
setwd("/Users/jz132/Desktop/hardac-xfer/MYC/p_value_list_fixed_analytical/enhancer/")
filenames <- mixedsort(Sys.glob("*.txt"))

n_enhancers <- nrow(data_enhancers_fantom)
enhancers_w_mutation <- which(data_enhancers_mutated$count > 0)

p_value_list <- rep(1, n_enhancers)
for(filename in filenames){
  new_p_value <- read_delim(filename, delim = " ", col_names = "p_value") %>% 
    pull(p_value)
  p_value_list <- pmin(p_value_list, new_p_value)
}

p_values_analytical <- p_value_list[enhancers_w_mutation]

# collect the simulation p_values result
setwd("/Users/jz132/Desktop/hardac-xfer/MYC/p_value_list_fixed_simulation/enhancer/")
filenames <- mixedsort(Sys.glob("*.txt"))

n_enhancers <- nrow(data_enhancers_fantom)
enhancers_w_mutation <- which(data_enhancers_mutated$count > 0)

p_value_list <- rep(1, n_enhancers)
for(filename in filenames){
  new_p_value <- read_delim(filename, delim = " ", col_names = "p_value") %>% 
    pull(p_value)
  p_value_list <- pmin(p_value_list, new_p_value)
}

p_values_simulation <- p_value_list[enhancers_w_mutation]

plot(p_values_analytical, p_values_simulation,
     main = "compare p-values between analytical and simulation analysis")

plot(-log10(p_values_analytical), -log10(p_values_simulation))
