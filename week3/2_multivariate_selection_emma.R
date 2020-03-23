# install.packages("tidyverse") # you only have to do this once on your personal computer
# install.packages("MCMCglmm") # you only have to do this once on your personal computer
library(tidyverse) # for data wrangling and plotting
library(MCMCglmm) # for linear regressions

# Example 1: Tameness and Tail Length are Associated. Ear Flop is independent.
# Generate phenotypes
set.seed(1)
tameness <- rnorm(1000, mean = 10, sd = 1)

set.seed(1)
tail_length <- rnorm(1000, mean = 10, sd = 1) + rnorm(1000, mean = 0.1, sd = 0.1)

set.seed(3)
ear_flop<- rnorm(1000, mean = 10, sd = 1)

# generate animal IDs
animal_id <- paste("animal",1:1000,sep="_")
phenotype_data <- cbind(tameness, tail_length, ear_flop)
row.names(phenotype_data) <- animal_id
head(phenotype_data)

# Generate a variance-covariance matrix for these phenotypes
# Note: this is a P matrix
# We are using it instead of a G matrix (matrix of additive genetic covariances and variances) for the sake of the example.
phenotype_matrix <- cov(phenotype_data)
phenotype_matrix

# Plot a heatmap of the associations
heatmap_colors <- colorRampPalette(c("blue","white", "red"))
heatmap(phenotype_matrix, col = heatmap_colors(100))

# Generate a selection gradient
head(phenotype_data)

# convert to data frame so we can add in other variables
phenotype_data_frame <- as.data.frame(phenotype_data)
phenotype_data_frame$animal_id <- row.names(phenotype_data_frame)
head(phenotype_data_frame)

# Add a vector for fitness. Associate fitness with tameness.
phenotype_data_frame$fitness <- ifelse(phenotype_data_frame$tameness>10,1,0)
summary(phenotype_data_frame)

# Look at the relationship between traits and fitness
plot(fitness~tameness,data=phenotype_data_frame)
plot(fitness~tail_length,data=phenotype_data_frame)
plot(fitness~ear_flop,data=phenotype_data_frame)

# Extract partial regression coefficients for each trait
set.seed(1)
selection_model <- MCMCglmm(fitness~tameness+tail_length+ear_flop,
                   random=~animal_id,
                   data=phenotype_data_frame)
summary(selection_model)

# Extract coefficients
selection_coeff_tameness <- mean(selection_model$Sol[,2])
selection_coeff_tail_length <- mean(selection_model$Sol[,3])
selection_coeff_ear_flop <- mean(selection_model$Sol[,4])

# create a selection gradient
selection_gradient <- t(cbind(selection_coeff_tameness,selection_coeff_tail_length, selection_coeff_ear_flop))
selection_gradient

# Calculate delta z
deltaz <- phenotype_matrix %*% selection_gradient
round(deltaz,3)

# What is the predicted change if selection acts ONLY on tameness?
# Set selection coefficient for other traits to ZERO
beta.onlyTameness <- selection_gradient
beta.onlyTameness[2:3] = 0
beta.onlyTameness
deltaz.onlyTameness = phenotype_matrix %*% beta.onlyTameness
round(deltaz.onlyTameness,3)
boxplot(t(deltaz.onlyTameness), main = "Expected Change if Selecting Only on Tameness", ylim=c(-1,1))
abline(h = 0, col = "red", lwd = 2)

# What is the predicted change if selection acts ONLY on tail length?
# Set selection coefficient for other traits to ZERO
beta.onlyTailLength <- selection_gradient
beta.onlyTailLength[-2] = 0
beta.onlyTailLength
deltaz.onlyTailLength = phenotype_matrix %*% beta.onlyTailLength
round(deltaz.onlyTailLength,3)
boxplot(t(deltaz.onlyTailLength), main = "Expected Change if Selecting Only on Tail Length", ylim=c(-1,1))
abline(h = 0, col = "red", lwd = 2)

# What is the predicted change if selection acts ONLY on Ear Flop?
# Set selection coefficient for other traits to ZERO
beta.onlyEarFlop <- selection_gradient
beta.onlyEarFlop[-3] = 0
beta.onlyEarFlop
deltaz.onlyEarFlop = phenotype_matrix %*% beta.onlyEarFlop
round(deltaz.onlyEarFlop,3)
boxplot(t(deltaz.onlyEarFlop), main = "Expected Change if Selecting Only on Ear Flop", ylim=c(-1,1))

abline(h = 0, col = "red", lwd = 2)

# Group 1 -------
# Work through an example where two traits have positive covariance, and one trait is under negative directional selection

phenotype_data_frame$fitness <- ifelse(phenotype_data_frame$tameness<10,1,0)
summary(phenotype_data_frame)

# Look at the relationship between traits and fitness
plot(fitness~tameness,data=phenotype_data_frame)
plot(fitness~tail_length,data=phenotype_data_frame)
plot(fitness~ear_flop,data=phenotype_data_frame)

# Extract partial regression coefficients for each trait
set.seed(1)
selection_model <- MCMCglmm(fitness~tameness+tail_length+ear_flop,
                            random=~animal_id,
                            data=phenotype_data_frame)
summary(selection_model)

# Extract coefficients
selection_coeff_tameness <- mean(selection_model$Sol[,2])
selection_coeff_tail_length <- mean(selection_model$Sol[,3])
selection_coeff_ear_flop <- mean(selection_model$Sol[,4])

# create a selection gradient
selection_gradient <- t(cbind(selection_coeff_tameness,selection_coeff_tail_length, selection_coeff_ear_flop))
selection_gradient

# Calculate delta z
deltaz <- phenotype_matrix %*% selection_gradient
round(deltaz,3)




# Group 2 --------
# Work through an example where two traits have positive covariance, and one trait is under positive directional selection and one trait is under negative directional selection

# Group 3 --------
# Work through an example where two traits have negative covariance, and one trait is under positive directional selection and one trait is under negative directional selection
