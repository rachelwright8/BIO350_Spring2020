library(tidyverse)

genos <- read.csv("nettie_family_data.csv")
head(genos)

# QUESTION 1: What range of IBS values would you expect between mother and daughter? ----
# plot IBS for Nettie and her mother
plot( jitter(2-abs( genos$nettie - genos$mother )), ylim=c(0,2), 
      ylab = "IBS", xlab = "SNP", main = "Mother", cex=0.1)

#QUESTION 2: What range of IBS values would you expect between full siblings? ----
# plot IBS for Nettie and her older brother
plot( jitter(2-abs( genos$nettie - genos$older_brother )), ylim=c(0,2), 
      ylab = "IBS", xlab = "SNP", main = "Nettie - Older Brother", cex=0.1)

#QUESTION 3: Who is most related? ----
# make functions to calculate IBS
IBSdiffs = function(x,y) round(mean(2-abs(x-y)),2) 
apply( genos, 2, IBSdiffs, y=genos$nettie)

# Plot IBS for Nettie and her younger brother
plot( jitter(2-abs( genos$nettie - genos$younger_brother )), ylim=c(0,2), 
      ylab = "IBS", xlab = "SNP", main = "Nettie - Younger Brother", cex=0.1)

#QUESTION 4: Nettie's older brother, father, grandfather, 
# and cousin share the same autosomal dominant disease. 
# Use these genotype data to eliminate regions of the chromosome 
# that are unlikely to carry the disease-causing locus.

# plot affected individuals
# linked region cannot be IBS0

#quartz()
#par(mfrow=c(2,2))
plot( jitter(2-abs( genos$grandfather - genos$older_brother )), ylim=c(0,2), 
      ylab = "IBS", xlab = "SNP", main = "Grandfather-Older Brother", cex=0.1)
plot( jitter(2-abs( genos$grandfather - genos$cousin_III2 )), ylim=c(0,2), 
      ylab = "IBS", xlab = "SNP", main = "Grandfather-Cousin", cex=0.1)
plot( jitter(2-abs( genos$father - genos$cousin_III2 )), ylim=c(0,2), 
      ylab = "IBS", xlab = "SNP", main = "Father-Cousin", cex=0.1)
plot( jitter(2-abs( genos$older_brother - genos$cousin_III2 )), ylim=c(0,2), 
      ylab = "IBS", xlab = "SNP", main = "Brother-Cousin", cex=0.1)
