load("FlowSorted.BloodExtended.EPIC.compTable.rda")#If EPIC
load("FlowSorted.BloodExtended.450klegacy.compTable")#If 450K
library(FlowSorted.Blood.EPIC)

#EPIC
Pred_EPIC <- projectCellType_CP(beta_matrix[rownames(FlowSorted.BloodExtended.EPIC.compTable),], 
                                         FlowSorted.BloodExtended.EPIC.compTable,lessThanOne =T)*100

#450K
Pred_450k <- projectCellType_CP(beta_matrix[rownames(FlowSorted.BloodExtended.450klegacy.compTable),], 
                           FlowSorted.BloodExtended.450klegacy.compTable,lessThanOne =T)*100

