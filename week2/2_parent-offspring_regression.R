install.packages("dplyr") # you only have to do this ONCE on your personal computer
library(dplyr) # but you must load the library every time you want to use it

# Generate Simulated Parental Phenotype Values ------
set.seed(1)
parents <- data.frame("pheno_dam" = rnorm(n = 1e3,mean = 10,sd = 1),
                      "pheno_sire" = rnorm(n = 1e3,mean = 11,sd = 1))
head(parents)

# Plot a Histogram of Parental Phenotypes
hist(parents$pheno_sire, col = "red", main = "Parent Phenotypes")
hist(parents$pheno_dam, col = "blue", add=T)

# Calculate Mid-Parent Values
parents <- parents %>% rowwise() %>% mutate(mid_parent_value = mean(c(pheno_dam,pheno_sire)))
head(parents)

# Make a vector of just mid-parent values to save ourselves some typing down the road
mpv <- parents$mid_parent_value

# Plot a Histogram of Mid-Parent Values
hist(parents$pheno_sire, col = "red", main = "Parent Phenotypes")
hist(parents$pheno_dam, col = "blue", add=T)
hist(mpv, col = "black", add=T)

# Generate Mean Offspring Values

# Simulate data where you know heritability will be extremely HIGH --------
pheno_offspring_high_heritability <- mpv

# Plot the Linear Regressions (High Heritability)
lm_high <- lm(mpv~pheno_offspring_high_heritability)
summary(lm_high)

plot(mpv, pheno_offspring_high_heritability, 
     col='purple', main = "Simulated High Heritability")

# The second coefficient from the object "coefficients" in "lm_high" is the slope
slope <- lm_high$coefficients[2] 
slope # this is our heritability

# Add the slope to the plot using the `abline` function
abline(lm_high, lwd=5)

# Another way to alculate heritability (same thing as slope above)
# Here we calculate heritability as the covariance between parents and offspring
# divided by the variance among parents
heritability_high <- cov(mpv,pheno_offspring_high_heritability)/var(mpv)
heritability_high

# Simulate data where you know heritability will be LOW ------
set.seed(1)
pheno_offspring_low_heritability <- sample(mpv)

# Low heritability
lm_low <- lm(mpv~pheno_offspring_low_heritability)
summary(lm_low)

plot(mpv, pheno_offspring_low_heritability, 
     col='green', main = "Simulated Low Heritability")

# The second coefficient from the object "coefficients" is the slope
slope <- lm_low$coefficients[2] 
slope

# Add the slope to the plot using the `abline` function
abline(lm_low, lwd=5)

# Another way to alculate heritability (same thing as slope above)
# Here we calculate heritability as the covariance between parents and offspring
# divided by the variance among parents
heritability_low <- cov(mpv,pheno_offspring_low_heritability)/var(mpv)
heritability_low

