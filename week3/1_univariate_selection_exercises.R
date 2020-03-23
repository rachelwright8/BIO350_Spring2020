#install.packages("tidyverse") # you only have to do this ONCE on your personal computer
library(tidyverse) # but you must load the library every time you want to use it

# Generate Simulated Parental Phenotype Values ------
set.seed(1)
parents <- data.frame("pheno_dam" = rnorm(1e3,10,1),
                      "pheno_sire" = rnorm(1e3,11,1))
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

# Add STRONG truncation selection ------------
selection_value <- 11.5 # only individuals with this value or higher will survive the selection event
hist(mpv)
abline(v = selection_value, col = "red", lwd = 3)

pheno_after_selection <- mpv[mpv>selection_value]
hist(pheno_after_selection, add = T, col = "red")

selection_differential <- mean(pheno_after_selection)-mean(mpv)
selection_differential # the difference in phenotypic means between all parents and parents that got a chance to reproduce

# apply the breeder's equation to predict the mean phenotypic change in the next generation
delta_z <- heritability_high*selection_differential
delta_z

# Calculate the mean in the next generation
offspring_after_selection <- mean(mpv)+delta_z

# Look at all of our values
mean_phenotypes <- data.frame("population" =  c("all_parents", "parents_after_selection", "offspring"), 
                              means = c(mean(mpv), mean(pheno_after_selection), offspring_after_selection),
                              sd = c(sd(mpv), sd(mpv), sd(mpv))) # Assume that variance remains constant (although it might not in reality!)

# change the order of the levels of "population" to make more sense when plotting
mean_phenotypes$population <- factor(mean_phenotypes$population,levels = c("all_parents", "parents_after_selection", "offspring"))

ggplot(mean_phenotypes, aes(x=population, y=means, group=1, col=population)) + 
  geom_line()+
  geom_pointrange(aes(ymin=means-sd, ymax=means+sd))+
  theme_bw()

# Plot this shift in mean phenotypic values
hist(mpv)
abline(v = mean(mpv), col = "black", lwd = 3)
abline(v = offspring_after_selection, col = "red", lwd = 3)

# Group 1-------------
# The trait has high heritability and selection is weak

head(mpv)
heritability <-  .8
selection_differential <- -.5  

delta_z <- heritability*selection_differential
delta_z

offspring_after_selection <- mean(mpv)+delta_z

pheno_after_selection <- selection_differential + mpv 

(head)pheno_after_selection

mean_phenotypes <- data.frame("population" =  c("all_parents", "parents_after_selection", "offspring"), 
                              means = c(mean(mpv), mean(pheno_after_selection), offspring_after_selection),
                              sd = c(sd(mpv), sd(mpv), sd(mpv))) # Assume that variance remains constant (although it might not in reality!)

mean_phenotypes$population <- factor(mean_phenotypes$population,levels = c("all_parents", "parents_after_selection", "offspring"))

ggplot(mean_phenotypes, aes(x=population, y=means, group=1, col=population)) + 
  geom_line()+
  geom_pointrange(aes(ymin=means-sd, ymax=means+sd))+
  theme_bw()


# Group 2 ------------
# The trait has low heritability and selection is strong

# Hint: start with the same parental values. Add your own values for heritability and selection.
head(mpv)
heritability <-  0.15
selection_differential <-  4
delta_z_Group2 <- heritability*selection_differential
offspring_after_selection_Group2 <- mean(mpv)+delta_z_Group2

pheno_after_selection_group2 <- mean(mpv)+selection_differential
mean_phenotypes <- data.frame("population" =  c("all_parents", "parents_after_selection", "offspring"), 
                              means = c(mean(mpv), pheno_after_selection_group2, offspring_after_selection_Group2),
                              sd = c(sd(mpv), sd(mpv), sd(mpv))) # Assume that variance remains constant (although it might not in reality!)
mean_phenotypes$population <- factor(mean_phenotypes$population,levels = c("all_parents", "parents_after_selection", "offspring"))

ggplot(mean_phenotypes, aes(x=population, y=means, group=1, col=population)) + 
  geom_line()+
  geom_pointrange(aes(ymin=means-sd, ymax=means+sd))+
  theme_bw()


# Group 3 -------------
# The trait has low heritability and selection is weak

head(mpv)
heritability <-  0.2
selection_differential <- 0.5 

# apply the breeder's equation to predict the mean phenotypic change in the next generation
delta_z <- heritability*selection_differential
delta_z

# Calculate the mean in the next generation
offspring_after_selection <- mean(mpv)+delta_z
pheno_after_selection <- selection_differential + mean(mpv)

# Look at all of our values
mean_phenotypes <- data.frame("population" =  c("all_parents", "parents_after_selection", "offspring"), 
                              means = c(mean(mpv), pheno_after_selection, offspring_after_selection),
                              sd = c(sd(mpv), sd(mpv), sd(mpv))) # Assume that variance remains constant (although it might not in reality!)

# change the order of the levels of "population" to make more sense when plotting
mean_phenotypes$population <- factor(mean_phenotypes$population,levels = c("all_parents", "parents_after_selection", "offspring"))

ggplot(mean_phenotypes, aes(x=population, y=means, group=1, col=population)) + 
  geom_line()+
  geom_pointrange(aes(ymin=means-sd, ymax=means+sd))+
  theme_bw()

# Plot this shift in mean phenotypic values
hist(mpv)
abline(v = mean(mpv), col = "black", lwd = 3)
abline(v = offspring_after_selection, col = "red", lwd = 3)


