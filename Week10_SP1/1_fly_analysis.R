setwd("~/Dropbox/BIO350/Week7_SP1_Flies/")

library(tidyverse)
library(MCMCglmm)
# install.packages("summarytools")
library(summarytools)

load("FlyData.Rdata")

# Table 1
view(dfSummary(phenotypes))

# Correlations Between Traits --------
plot(PostMatingEgg ~ PreMatingEgg, data = phenotypes)
lm_egg_production <- lm(PostMatingEgg ~ PreMatingEgg, data = phenotypes)
abline(lm_egg_production, col = "red")
summary(lm_egg_production)
cor(phenotypes$PostMatingEgg, phenotypes$PreMatingEgg, 
    method = "pearson", use = "complete.obs")

plot(Duration ~ PreMatingEgg, data = phenotypes)
lm_mating_duration <- lm(Duration ~ PreMatingEgg, data = phenotypes)
abline(lm_mating_duration, col = "red")
summary(lm_mating_duration)

plot(PostMatingEgg ~ Duration, data = phenotypes)
lm_egg_mating <- lm(PostMatingEgg ~ Duration, data = phenotypes)
abline(lm_egg_mating, col = "red")
summary(lm_egg_mating)

# Group 1: Egg Production Post-Mating --------
# Partition variance using MCMCglmm
head(pedigree)
head(phenotypes)

# set the "prior" for the model
prior <- list(G = list(G1 = list(V = 1, n = 1), 
                        G2 = list(V = 1, n = 1)), 
               R = list(V = 10, n = 1))

# Run the MCMC model
set.seed(1)
model_EPpost <- MCMCglmm(PostMatingEgg ~ Block, random = ~animal+MaleID, 
                           pedigree = pedigree, data = phenotypes, 
                           prior = prior)

# Additive genetic variance explained by female ID (animal) and male ID
summary(model_EPpost)$Gcovariances

Va_female <- summary(model_EPpost)$Gcovariances[1]
Va_female
Va_male <- summary(model_EPpost)$Gcovariances[2]
Va_male

# Residual variance
summary(model_EPpost)$Rcovariances

# Total phenotypic variance
Vtotal <- summary(model_EPpost)$Gcovariances[1]+
  summary(model_EPpost)$Gcovariances[2]+
  summary(model_EPpost)$Rcovariances[1]
Vtotal

# Heritability

# Female component
h2_female <- Va_female/Vtotal
h2_female

# Male component
h2_male <- Va_male/Vtotal
h2_male

# Total heritability
tau <- h2_female+h2_male
tau

# Group 2: Copulation Duration --------
# set the "prior" for the model
prior <- list(G = list(G1 = list(V = 1, n = 1), 
                        G2 = list(V = 1, n = 1)), 
               R = list(V = 10, n = 1))

# Run the MCMC model
set.seed(1)
model_duration <- MCMCglmm(Duration ~ Block, random = ~animal+MaleID, 
                           pedigree = pedigree, data = phenotypes, family="gaussian",
                           prior = prior)

# Additive genetic variance explained by female ID (animal) and male ID
summary(model_duration)$Gcovariances

Va_female <- summary(model_duration)$Gcovariances[1]
Va_female
Va_male <- summary(model_duration)$Gcovariances[2]
Va_male

# Residual variance
summary(model_duration)$Rcovariances

# Total phenotypic variance
Vtotal <- summary(model_duration)$Gcovariances[1]+
  summary(model_duration)$Gcovariances[2]+
  summary(model_duration)$Rcovariances[1]
Vtotal

# Heritability

# Female component
h2_female <- Va_female/Vtotal
h2_female

# Male component
h2_male <- Va_male/Vtotal
h2_male

# Total heritability
tau <- h2_female+h2_male
tau

# Group 3: Egg Production Pre-Mating --------
# Partition variance using MCMCglmm
head(pedigree)
head(phenotypes)

# set the "prior" for the model
prior <- list(G = list(G1 = list(V = 1, n = 1), 
                       G2 = list(V = 1, n = 1)), 
              R = list(V = 10, n = 1))

# Run the MCMC model
set.seed(1)
model_EPpre <- MCMCglmm(PreMatingEgg ~ Block, random = ~animal+MaleID, 
                         pedigree = pedigree, data = phenotypes,family="poisson", 
                         prior = prior)

# Additive genetic variance explained by female ID (animal) and male ID
summary(model_EPpre)$Gcovariances

Va_female <- summary(model_EPpre)$Gcovariances[1]
Va_female
Va_male <- summary(model_EPpre)$Gcovariances[2]
Va_male

# Residual variance
summary(model_EPpre)$Rcovariances

# Total phenotypic variance
Vtotal <- summary(model_EPpre)$Gcovariances[1]+
  summary(model_EPpre)$Gcovariances[2]+
  summary(model_EPpre)$Rcovariances[1]
Vtotal

# Heritability

# Female component
h2_female <- Va_female/Vtotal
h2_female

# Evolve the poluation --------

# Phenotypic co-variance
head(phenotypes)
phenotype_data <- phenotypes %>% select(PreMatingEgg, Duration, PostMatingEgg) %>% filter(complete.cases(.))
head(phenotype_data)
phenotype_matrix <- cov(phenotype_data)
phenotype_matrix

# create a selection gradient
selection_gradient <- c(1,-1,0)
selection_gradient

# Calculate delta z
deltaz <- phenotype_matrix %*% selection_gradient
round(deltaz,3)


