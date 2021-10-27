library(gtools)
setwd("/Users/jz132/Desktop/hardac-xfer/p_value_list_fixed/enhancer")
filenames <- mixedsort(Sys.glob("*.txt"))

temp_filename <- filenames[1]
n_enhancers <- nrow(read.table(temp_filename))

p_value_list_fixed <- numeric(n_enhancers)
for(filename in filenames){
  new_p_value <- read.table(filename)[,1]
  p_value_list_fixed <- p_value_list_fixed + new_p_value
}

p_values_fixed <- p_value_list_fixed/length(filenames)
p_values_fixed <- (p_values_fixed*10^6 + 1)/(10^6 + 1) # correct the monte carlo p-values

# p_values_fixed_non_one <- p_values_fixed
p_values_fixed_non_one <- p_values_fixed[p_values_fixed < 1]
# y <- qunif(ppoints(length(p_values_fixed_non_one)))
# qqplot(y, p_values_fixed_non_one,
#        main = "enhancer mutation p-values fixed number of mutations")
# abline(a = 0, b= 1)

y <- qunif(ppoints(length(p_values_fixed_non_one)))
qqplot(-log(y, base = 10), -log(p_values_fixed_non_one, base = 10),
       main = "enhancer mutation p-values fixed number of mutations")
abline(a = 0, b= 1)

hist(p_values_fixed)
hist(p_values_fixed_non_one)
