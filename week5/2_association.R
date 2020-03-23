install.packages("rrBLUP")
library(rrBLUP)
# This is all based on the rrBLUP vignette
# Endelman JB (2011). “Ridge regression and other kernels for genomic selection with R package rrBLUP.” Plant Genome, 4, 250-255.

markers <- read.csv("call_method_32.b",skip=1,header=T,as.is=T,check.names=FALSE)
markers[1:5,1:5]
dim(markers)

# create a function to convert ATGC format to {-1,0,1,NA}
convert.snp  <- function(x) {
    alleles <- na.omit(unique(x))
    y <- rep(NA,length(x))
    y[which(x==alleles[1])] <- -1
    y[which(x==alleles[2])] <- 1
    return(y)
}

M <- apply(markers[,-(1:2)],1,convert.snp)
M[1:5,1:5]
dim(M)
# number of lines number of markers
# 199             216130

# create genotype identifiers (gid)
gid <- colnames(markers)[-(1:2)]
rownames(M) <- gid
n <- nrow(M) # number of lines
m <- ncol(M) # number of markers
head(gid)

# create a marker-based relationship matrix
# genetic varianace and covariance among lines
A <- A.mat(M)
A[1:5,1:5]

# The GWAS function we will use expects the genotype data to have columns for "marker", "chromosome", "position" and then all of the genetic lines
geno <- cbind(1:m,markers[,1:2],t(M))
colnames(geno) <- c("marker","chrom","pos",gid)
geno[1:10,1:10]

# Load in the phenotype data
pheno0 <- read.table("phenotype_published_raw.tsv",header=T,as.is=T,check.names=FALSE,sep="\t")
pheno0[1:5,1:5]

# Extract three phenotypes
# FT_LD	= flowering	time under long	days
# Dormancy	=	seed dormancy
# avrRpm	=	disease	resistance
pheno <- pheno0[,c(1,3,10,35)]
pheno[1:5,1:4]
# Make the names prettier
colnames(pheno) <- c("ecoid","FloweringTime","Dormancy","Resistance")
head(pheno)

group1 <- pheno[c(1,2)] # just flowering time
group2 <- pheno[c(1,3)] # just dormancy
group3 <- pheno[c(1,4)] # just disease resistance

# Calculate GWAS
# This will take awhile
gwas_flowering <- GWAS(pheno=group1, # your phenotype data
                  geno=geno, # your genotype data
                  K=A) # a kinship matrix

gwas_dormancy <- GWAS(pheno=group2, # your phenotype data
                       geno=geno, # your genotype data
                       K=A) # a kinship matrix

gwas_resistance <- GWAS(pheno=group3, # your phenotype data
                       geno=geno, # your genotype data
                       K=A) # a kinship matrix

# Reduce the sample size
reduced_n <- 50
reduced_M <- M[1:reduced_n,]
dim(reduced_M)
reduced_M[1:5,1:5]
reduced_A <- A.mat(reduced_M)
reduced_A[1:5,1:5]

# The GWAS function we will use expects the genotype data to have columns for "marker", "chromosome", "position" and then all of the genetic lines
reduced_geno <- cbind(1:m,markers[,1:2],t(reduced_M))
colnames(reduced_geno) <- c("marker","chrom","pos",gid[1:reduced_n])
reduced_geno[1:10,1:10]
dim(reduced_geno)

reduced_pheno <- pheno[pheno$ecoid %in% names(reduced_geno),]
dim(pheno)
dim(reduced_pheno)

reduced_group1 <- reduced_pheno[c(1,2)] # just flowering time
reduced_group2 <- reduced_pheno[c(1,3)] # just dormancy
reduced_group3 <- reduced_pheno[c(1,4)] # just disease resistance

# Calculate GWAS
# This will take awhile
gwas_flowering <- GWAS(pheno=reduced_group1, # your phenotype data
                       geno=reduced_geno, # your genotype data
                       K=reduced_A) # a kinship matrix

gwas_dormancy <- GWAS(pheno=reduced_group2, # your phenotype data
                      geno=reduced_geno, # your genotype data
                      K=reduced_A) # a kinship matrix

gwas_resistance <- GWAS(pheno=reduced_group3, # your phenotype data
                        geno=reduced_geno, # your genotype data
                        K=reduced_A) # a kinship matrix
