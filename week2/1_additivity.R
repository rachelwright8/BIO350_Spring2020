#install.packages("ggplot2") # you only have to do this ONCE on your personal computer
library(ggplot2)

# Set allele frequencies, phenotypic values, and dominance ---------

# allele frequency for `A`
p <- 0.15

# allele frequency for `a`
q <- 1-p
q

# phenotypic value for `AA`
a <- 1

# phenotypic value for `aa`
-a

# phenotypic value for `Aa`
d <- 1

# dominance
d/a

# Calculate the difference between the mean phenotype for 
# each genotype and the population mean
# Some tedious algebra comes out to...
alpha <- a+((q-p)*d)
alpha

# Play with the numbers above to see if this makes sense.

# Calculate genotypic variance ---------

# total genotypic variance is
Var_genotypic <- 2*p*q*alpha^2 + (2*p*q*d)^2

# then additive genetic variance is
Var_additive <- 2*p*q*alpha^2

# and dominance variance is
Var_dominance <- (2*p*q*d)^2

# Create a data frame of the results

phenotypic_variance <- data.frame("types"=c("total", "additive", "dominance"),
                                  "values"=c(Var_genotypic, Var_additive, Var_dominance))
phenotypic_variance 

# Plot Results

ggplot(data=phenotypic_variance, aes(x = types, y = values)) +
  geom_bar(stat="identity") + 
  theme_bw()
