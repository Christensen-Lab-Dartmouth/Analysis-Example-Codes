library(sesame)
library(sesameData)
workdir<-"." #Direct to idats folder
filtered_betas <- openSesame(workdir, func = getBetas, prep = "QCDPB")#Change QCDPB to TQCDPB for mouse
unfiltered_betas<-openSesame(workdir, func = getBetas, prep="QCDPB", mask = FALSE)#Change QCDPB to TQCDPB for mouse
                               