BiocManager::install("FlowSorted.Blood.EPIC")
library(FlowSorted.Blood.EPIC)
library(minfi)

RGset <- read.metharray.exp("Directory to idats",force = TRUE)

MSet_noob <- preprocessNoob(RGset)
beta_matrix<-getBeta(MSet_noob)


#EPIC
Pred_EPIC <- projectCellType_CP(beta_matrix[rownames(IDOLOptimizedCpGs.compTable),], 
                                IDOLOptimizedCpGs.compTable,lessThanOne =T)*100

#450K
Pred_450K <- projectCellType_CP(beta_matrix[rownames(IDOLOptimizedCpGs450klegacy.compTable),], 
                                IDOLOptimizedCpGs450klegacy.compTable,lessThanOne =T)*100