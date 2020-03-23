# Estimate broad-sense heritability using animal models

# Generate a vector of genotypes (5 clones per genotype)
gtypes <- rep(letters,each=5)
gtypes

# Generate a vector of phenotypes for each individual
# Start with a scenario where heritability is HIGH ------
# i.e., clones are very similar to one another
set.seed(1)
ptypes <- rep(rnorm(26),each=5) + (rnorm(5*26))
head(ptypes)
# What is the total observed variance?
var(ptypes)

# Partition the variance due to genetics and residiual effects using ANOVA
av <- anova(lm(ptypes~gtypes))
av

# Extract the residual variance
res_variance <- av[[3]][2] # the value in the third column, second row

# Subtract the residual variance from variance attributed to genotype
# divide by the number of individuals per genotype
gen_variance <- (av[[3]][1]-res_variance)/5

# Calculate broad sense heritability
broad_sense <- gen_variance/(gen_variance+res_variance)
broad_sense

# Use a package to do the same thing!
install.packages("heritability") # you only have to do this ONCE on your personal computer
library(heritability)

repeatability(data.vector = ptypes, geno.vector = gtypes)

# Do the values match?

# Now try a scenario where heritability is LOW -------
# i.e., clones are NOT very similar to one another
set.seed(1)
ptypes <- rep(rnorm(26),each=5) + (rnorm(5*26))*3

# What is the total observed variance?
var(ptypes)

# Partition the variance due to genetics and residiual effects using ANOVA
av <- anova(lm(ptypes~gtypes))
av

# Extract the residual variance
res_variance <- av[[3]][2]

# Subtract the residual variance from variance attributed to genotype
# divide by the number of individuals per genotype
gen_variance <- (av[[3]][1]-res_variance)/5

# Calculate broad sense heritability
broad_sense <- gen_variance/(gen_variance+res_variance)
broad_sense

# Using the heritability package to calculate LOW heritability
repeatability(data.vector = ptypes, geno.vector = gtypes)

# Real Data from the Paper -----------

# Load the data from the Jury, Delano, and Toonen paper
coral_weights <- read.table("coral_data.txt", header=T)

# Make sure it looks OK
head(coral_weights)
str(coral_weights)

# Notice that there are multiple covariates (possible sources of variation) that need to be considered by the model

repeatability(data.vector = coral_weights$growthRate, 
              geno.vector = coral_weights$colony,
              covariates.frame = data.frame("treatment" = coral_weights$pH,
                                            "collection" = coral_weights$collectionSite), 
              line.repeatability = F)
