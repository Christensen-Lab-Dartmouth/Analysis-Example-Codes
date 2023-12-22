library(limma)
library(minfi)
library(dplyr)


# Assume you have a methylation beta matrix and Pheno matrix which contains phenotype variables for all samples
# Run epigenome-wide association study (EWAS) to identify CpGs that are differentially methylated CpGs between case and control

CaseStatus <- factor(Pheno$CaseStatus)#Assume CaseStatus is a binary variable with Yes and No in Pheno

#Assume you want to control for age and sex as confounders in your model
Age<-Pheno$Age
Sex<-factor(Pheno$Sex)


# use the above to create a design matrix
design <- model.matrix(~0+CaseStatus+Age+Sex)#adjusted for age,sex,tumor,CD8mem,DC,Epithelial
colnames(design) <- c(levels(CaseStatus),"Age",levels(Sex)[-1])

library(limma)
fit1 <- lmFit(EWAS_Beta, design)

contMatrix1 <- makeContrasts(Yes-No,levels = design)

fit2 <- contrasts.fit(fit1, contMatrix1)
fit2 <- eBayes(fit2)


# look at the numbers of differentially methylated CpGs at FDR < 0.05
summary(decideTests(fit2,adjust.method = "fdr",p.value = 0.05))

