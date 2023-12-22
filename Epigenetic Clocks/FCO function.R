#####################################################################################
#####################################################################################
# Title: "Tracing human stem cell lineage during development using DNA methylation" 
# Authors: Salas LA, Wiencke JK, Koestler DC, Christensen BC, Zhang Z and Kelsey KT
# 2018
#####################################################################################

#####################################################################################
# FUNCTION:  PredictFetalSignature 
#    This function predicts the fraction of cells carrying the fetal/ESC signature in 
#    blood-derived DNA methylation data
#
# ARGUMENTS:
#    Y:              Data frame (J x N) of methylation beta values for target 
#                    data set, i.e., blood-based DNAm data
#
#
# RETURNS:   A N x 2 matrix whose columns represent the predicted fraction of cells
#            carrying the Adult and Fetal/ESC signatures, respectively.
#
# USAGE:     In order to use the PredictFetalSignature function, the user first needs 
#            needs to first load "PredictFetalSignature.RData", which contains the 
#            27 CpG signature needed for deconvolution.  Further details on the 
#            application of the PredictFetalSignature function are given in the
#            Supplementary Materials
#####################################################################################
library(tibble)
library(openxlsx)
load("PredictFetalSignature.RData") #load FCO data and function
df <- read.table("your_data.csv", header = TRUE, sep = ",", stringsAsFactors = F, na.strings = "NA") #read your csv methylation beta data
betas <- data.frame(df[,-1], row.names=df[,1]) #clean data
FCO_Prediction<-PredictFetalSignature(betas)#run FCO
FCO_Prediction <- rownames_to_column(as.data.frame(FCO_Prediction))


