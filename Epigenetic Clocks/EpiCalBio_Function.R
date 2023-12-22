library(sesame)#verison 1.14.2
library(minfi)
library(ENmix)
library(FlowSorted.Blood.EPIC)
library(tibble)
library(dplyr)
EpiBioCal <- function(workdir) {
  sdfs = openSesame(workdir, func = NULL)
  data_pOOBHA<-as.data.frame(matrix(NA,nrow = nrow(sdfs[[1]]), ncol = length(sdfs)))
  rownames(data_pOOBHA)<-sdfs[[1]]$Probe_ID
  colnames(data_pOOBHA)<-names(sdfs)
  for (i in 1:length(sdfs)) {
    data_pOOBHA[,i]<-pOOBAH(sdfs[[i]], return.pval=TRUE)
  }
  #--------------------------------------------------------------
  data_sex<-as.data.frame(matrix(NA,nrow = length(sdfs) , ncol = 1))
  rownames(data_sex)<-names(sdfs)
  colnames(data_sex)<-"InferredSex"
  for (i in 1:length(sdfs)) {
    data_sex[i,1]<-inferSex(sdfs[[i]])
  }  
  data_sex<-rownames_to_column(data_sex)
  colnames(data_sex)[1]<-"Sample_ID"
  #--------------------------------------------------
  data_ethnicity<-as.data.frame(matrix(NA,nrow = length(sdfs) , ncol = 1))
  rownames(data_ethnicity)<-names(sdfs)
  colnames(data_ethnicity)<-"InferredEthnicity"
  for (i in 1:length(sdfs)) {
    data_ethnicity[i,1]<-inferEthnicity(sdfs[[i]])
  }  
  data_ethnicity<-rownames_to_column(data_ethnicity)
  colnames(data_ethnicity)[1]<-"Sample_ID"
  #-----------------------------------------------------------------
  data_RGset <- read.metharray.exp(workdir,force = TRUE)
  data_betas <- preprocessNoob(data_RGset)
  data_betas<-getBeta(data_betas)
  #load clocks
  load("clock_list.RDATA")
  output<-as.data.frame(matrix(NA,length(sdfs),7))
  colnames(output)<-c("Dunedin","Horvath","Hannum","PhenoAge","ZhangEN","EpiTOC2_TNSC","FCO")
  rownames(output)<-names(sdfs)
  #Dunedin----------------------------------------
  data_pOOBHA_Dunedin<-data_pOOBHA[clock_list[clock_list$Biomarker == "Dunedin","CpG"],]
  data_betas_Dunedin<-data_betas[,colnames(data_betas[,colSums(data_pOOBHA_Dunedin>0.05)<dim(data_pOOBHA_Dunedin)[1]*0.1])]
  Dunedin_Age<-methyAge(data_betas_Dunedin, normalize = F,nCores = 1)
  rownames(Dunedin_Age)<-Dunedin_Age$SampleID
  for (i in rownames(Dunedin_Age)) {
    output[i,"Dunedin"]<-Dunedin_Age[i,"PACE"]
  }
  #Horvath-------------------------------------------
  diff<-setdiff(clock_list[clock_list$Biomarker == "Horvath","CpG"],rownames(data_betas))
  data_pOOBHA_Horvath<-data_pOOBHA[rownames(data_pOOBHA)%in%clock_list[clock_list$Biomarker == "Horvath","CpG"],]
  data_betas_Horvath<-data_betas[,colnames(data_betas[,(colSums(data_pOOBHA_Horvath>0.05)+length(diff))<dim(data_pOOBHA_Horvath)[1]*0.1])]
  Horvath_Age<-methyAge(data_betas_Horvath, normalize = F,nCores = 1)
  rownames(Horvath_Age)<-Horvath_Age$SampleID
  for (i in rownames(Horvath_Age)) {
    output[i,"Horvath"]<-Horvath_Age[i,2]
  }
  #Hannum-------------------------------------------
  diff<-setdiff(clock_list[clock_list$Biomarker == "Hannum","CpG"],rownames(data_betas))
  data_pOOBHA_Hannum<-data_pOOBHA[rownames(data_pOOBHA)%in%clock_list[clock_list$Biomarker == "Hannum","CpG"],]
  data_betas_Hannum<-data_betas[,colnames(data_betas[,(colSums(data_pOOBHA_Horvath>0.05)+length(diff))<dim(data_pOOBHA_Hannum)[1]*0.1])]
  Hannum_Age<-methyAge(data_betas_Hannum, normalize = F,nCores = 1)
  rownames(Hannum_Age)<-Hannum_Age$SampleID
  for (i in rownames(Hannum_Age)) {
    output[i,"Hannum"]<-Hannum_Age[i,3]
  }
  #PhenoAge-------------------------------------------
  data_pOOBHA_PhenoAge<-data_pOOBHA[clock_list[clock_list$Biomarker == "PhenoAge","CpG"],]
  data_betas_PhenoAge<-data_betas[,colnames(data_betas[,colSums(data_pOOBHA_PhenoAge>0.05)<dim(data_pOOBHA_PhenoAge)[1]*0.1])]
  PhenoAge_Age<-methyAge(data_betas_PhenoAge, normalize = F,nCores = 1)
  rownames(PhenoAge_Age)<-PhenoAge_Age$SampleID
  for (i in rownames(PhenoAge_Age)) {
    output[i,"PhenoAge"]<-PhenoAge_Age[i,4]
  }
  #ZhangEN-------------------------------------------
  data_pOOBHA_ZhangEN<-data_pOOBHA[clock_list[clock_list$Biomarker == "Zhang_EN","CpG"],]
  data_betas_ZhangEN<-data_betas[,colnames(data_betas[,colSums(data_pOOBHA_ZhangEN>0.05)<dim(data_pOOBHA_ZhangEN)[1]*0.1])]
  
  read.table("en.coef",stringsAsFactor=F,header=T)->encoef
  en_int<-encoef[1,2]
  encoef<-encoef[-1,]
  rownames(encoef)<-encoef$probe
  dataNona<-t(data_betas_ZhangEN)
  dataNona.norm<- apply(dataNona,1,scale)  
  rownames(dataNona.norm)<-colnames(dataNona)
  encomm<- intersect(rownames(encoef),rownames(dataNona.norm))
  encoef<-encoef[encomm,]
  encoef$coef%*%dataNona.norm[encomm,]+en_int->enpred
  enpred<-t(enpred)
  colnames(enpred)<-"ZhangEN"
  for (i in rownames(enpred)) {
    output[i,"ZhangEN"]<-enpred[i,1]
  }
  #EpiTOC2------------------------------------------------------------
  load("dataETOC2.Rd")
  EpiTOC2 <- function(data.m,ages.v=NULL){
    # load("/dartfs/rc/lab/S/SalasLab/Rotations/Irma_Vlasac/Processed_data/GSE99553/EpiTOC2/dataETOC2.Rd"); ## this loads the CpG information
    cpgETOC.v <- dataETOC2.l[[2]];
    estETOC2.m <- dataETOC2.l[[1]];
    soloCpG.v <- dataETOC2.l[[3]];
    ### do epiTOC
    common.v <- intersect(rownames(data.m),cpgETOC.v);
    print(paste("Number of represented epiTOC CpGs (max=385)=",length(common.v),sep=""));
    map.idx <- match(common.v,rownames(data.m));
    pcgtAge.v <- colMeans(data.m[map.idx,],na.rm=TRUE);
    ### do EpiTOC2
    map.idx <- match(rownames(estETOC2.m),rownames(data.m));
    rep.idx <- which(is.na(map.idx)==FALSE);
    print(paste("Number of represented EpiTOC2 CpGs (max=163)=",length(rep.idx),sep=""))
    tmp.m <- data.m[map.idx[rep.idx],];
    TNSC.v <- 2*colMeans(diag(1/(estETOC2.m[rep.idx,1]*(1-estETOC2.m[rep.idx,2]))) %*% (tmp.m - estETOC2.m[rep.idx,2]),na.rm=TRUE);
    TNSC2.v <- 2*colMeans(diag(1/estETOC2.m[rep.idx,1]) %*% tmp.m,na.rm=TRUE);
    ### do HypoClock
    common.v <- intersect(rownames(data.m),soloCpG.v);
    print(paste("Number of represented solo-WCGWs (max=678)=",length(common.v),sep=""));
    map.idx <- match(common.v,rownames(data.m));
    hypoSC.v <- colMeans(data.m[map.idx,],na.rm=TRUE);
    
    estIR.v <- NULL; estIR2.v <- NULL;
    estIR <- NULL;  estIR2 <- NULL;
    if(!is.null(ages.v)){
      estIR.v <- TNSC.v/ages.v;
      estIR <- median(estIR.v,na.rm=TRUE);
      estIR2.v <- TNSC2.v/ages.v;
      estIR2 <- median(estIR2.v,na.rm=TRUE);
    }
    
    
    return(list(tnsc=TNSC.v,tnsc2=TNSC2.v,irS=estIR.v,irS2=estIR2.v,irT=estIR,irT2=estIR2,pcgtAge=pcgtAge.v,hypoSC=hypoSC.v));
  }
  diff<-setdiff(clock_list[clock_list$Biomarker == "EpiTOC2","CpG"],rownames(data_betas))
  data_pOOBHA_EpiTOC2<-data_pOOBHA[rownames(data_pOOBHA)%in%clock_list[clock_list$Biomarker == "EpiTOC2","CpG"],]
  data_betas_EpiTOC2<-data_betas[,colnames(data_betas[,(colSums(data_pOOBHA_EpiTOC2>0.05)+length(diff))<dim(data_pOOBHA_EpiTOC2)[1]*0.1])]
  epitoc_output <- EpiTOC2(data.m = data_betas_EpiTOC2)
  epitoc_tnsc<-as.data.frame(epitoc_output$tnsc2)
  colnames(epitoc_tnsc)<-"EpiTOC2"
  for (i in rownames(epitoc_tnsc)) {
    output[i,"EpiTOC2_TNSC"]<-epitoc_tnsc[i,1]
  }
  #FCO---------------------------------------------------------------
  load("PredictFetalSignature.RData")
  PredictFetalSignature<-function(Y) {
    
    projectWBCnew = function(Y, coefWBC, contrastWBC=NULL, nonnegative=TRUE){ 
      if(is.null(contrastWBC)) Xmat = coefWBC
      else Xmat = coefWBC %*% t(contrastWBC) 
      
      nCol = dim(Xmat)[2]
      nSubj = dim(Y)[2]
      
      mixCoef = matrix(0, nSubj, nCol)
      rownames(mixCoef) = colnames(Y)
      colnames(mixCoef) = colnames(Xmat)
      
      if(nonnegative){
        library(quadprog)
        
        Amat = cbind(rep(-1,nCol), diag(nCol))
        b0vec = c(-1,rep(0,nCol))
        
        for(i in 1:nSubj){
          obs = which(!is.na(Y[,i])) 
          Dmat = t(Xmat[obs,])%*%Xmat[obs,]
          mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,])%*%Y[obs,i], Amat, b0vec, meq = 0)$sol
        }
      }
      
      else {
        for(i in 1:nSubj){
          obs = which(!is.na(Y[,i])) 
          Dmat = t(Xmat[obs,])%*%Xmat[obs,]
          mixCoef[i,] = solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
        }
      }
      return(mixCoef)
    }
    
    invarCpGs = rownames(meanMethSignature)
    int = intersect(rownames(Y), invarCpGs)
    overlap = length(int)
    print(paste("Of the 27 invariant CpGs, ", overlap, " out of 27 were contained in the supplied data set", 
                sep = ""))
    targetData = Y[int,]
    projData = meanMethSignature[int,]
    pred = projectWBCnew(targetData,  projData)
    pred0 = ifelse(sign(pred) == -1, 0, pred)
    pred1 = t(apply(pred0, 1, function(w) w/sum(w)))
    pred1*100
  }
  diff<-setdiff(clock_list[clock_list$Biomarker == "FCO","CpG"],rownames(data_betas))
  data_pOOBHA_FCO<-data_pOOBHA[rownames(data_pOOBHA)%in%clock_list[clock_list$Biomarker == "FCO","CpG"],]
  data_betas_FCO<-data_betas[,colnames(data_betas[,(colSums(data_pOOBHA_FCO>0.05)+length(diff))<dim(data_pOOBHA_FCO)[1]*0.1])]
  FCO_Prediction<-PredictFetalSignature(data_betas_FCO)
  for (i in rownames(FCO_Prediction)) {
    output[i,"FCO"]<-FCO_Prediction[i,2]
  }
  if (nrow(sdfs[[1]])<800000){
    load("FlowSorted.BloodExtended.450klegacy.compTable.rda")
    
    data_pOOBHA_ED<-data_pOOBHA[rownames(FlowSorted.BloodExtended.450klegacy.compTable),]
    data_betas_ED<-data_betas[,colnames(data_betas[,colSums(data_pOOBHA_ED>0.05)<dim(data_pOOBHA_ED)[1]*0.1])]
    
    Pred_450k_OG <- projectCellType_CP(data_betas_ED[rownames(FlowSorted.BloodExtended.450klegacy.compTable),], 
                                       FlowSorted.BloodExtended.450klegacy.compTable,lessThanOne =T)*100
    
    Pred_450k_scaled<- Pred_450k_OG
    for (i in 1:nrow(Pred_450k_scaled)) {
      Pred_450k_scaled[i,]<-(Pred_450k_scaled[i,]/rowSums(Pred_450k_scaled)[i])*100
    }
    colnames(Pred_450k_scaled)<-paste0(colnames(Pred_450k_scaled),"2")
    Pred_450k<-as.data.frame(cbind(Pred_450k_OG,Pred_450k_scaled))
    
    Pred_450k<-rownames_to_column(Pred_450k)
    output<-rownames_to_column(output)
    output<-left_join(output,Pred_450k,by = "rowname")
    colnames(output)[1]<-"Sample_ID"
  }else{
    load("FlowSorted.BloodExtended.EPIC.compTable.rda")
    
    data_pOOBHA_ED<-data_pOOBHA[rownames(FlowSorted.BloodExtended.EPIC.compTable),]
    data_betas_ED<-data_betas[,colnames(data_betas[,colSums(data_pOOBHA_ED>0.05)<dim(data_pOOBHA_ED)[1]*0.1])]
    
    Pred_EPIC_OG <- projectCellType_CP(data_betas_ED[rownames(FlowSorted.BloodExtended.EPIC.compTable),], 
                                       FlowSorted.BloodExtended.EPIC.compTable,lessThanOne =T)*100
    
    Pred_EPIC_scaled<- Pred_EPIC_OG
    for (i in 1:nrow(Pred_EPIC_scaled)) {
      Pred_EPIC_scaled[i,]<-(Pred_EPIC_scaled[i,]/rowSums(Pred_EPIC_scaled)[i])*100
    }
    colnames(Pred_EPIC_scaled)<-paste0(colnames(Pred_EPIC_scaled),"2")
    Pred_EPIC<-as.data.frame(cbind(Pred_EPIC_OG,Pred_EPIC_scaled))
    
    Pred_EPIC<-rownames_to_column(Pred_EPIC)
    output<-rownames_to_column(output)
    output<-left_join(output,Pred_EPIC,by = "rowname")    
    colnames(output)[1]<-"Sample_ID"
  }
  output<-left_join(output,data_sex,by = "Sample_ID")
  output<-left_join(output,data_ethnicity,by = "Sample_ID")
  return(output)
}
workdir<-"test_450k"
test_data<-EpiBioCal(workdir = workdir)
workdir<-"test_epic"
test_data<-EpiBioCal(workdir = workdir)
workdir<-"GSE174555_RAW"
test_data<-EpiBioCal(workdir = workdir)
