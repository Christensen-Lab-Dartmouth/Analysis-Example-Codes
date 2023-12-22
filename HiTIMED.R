library(sesame)
library(sesameData)
library(dplyr)
library(minfi)
library(FlowSorted.Blood.EPIC)
library(HiTIMED)


Amir_betas<-openSesame("IDATS/", func = getBetas, prep="QCDPB", mask = FALSE)

HiTIMED_result<-HiTIMED_deconvolution(Amir_betas,"UCEC",6,"tumor")

