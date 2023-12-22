
load("dataETOC2.Rd")
epiTOC2 <- function(data.m,ages.v=NULL){
 # load("/dartfs/rc/lab/S/SalasLab/Rotations/Irma_Vlasac/Processed_data/GSE99553/EpiTOC2/dataETOC2.Rd"); ## this loads the CpG information
  cpgETOC.v <- dataETOC2.l[[2]];
  estETOC2.m <- dataETOC2.l[[1]];
  soloCpG.v <- dataETOC2.l[[3]];
  ### do epiTOC
  common.v <- intersect(rownames(data.m),cpgETOC.v);
  print(paste("Number of represented epiTOC CpGs (max=385)=",length(common.v),sep=""));
  map.idx <- match(common.v,rownames(data.m));
  pcgtAge.v <- colMeans(data.m[map.idx,],na.rm=TRUE);
  ### do epiTOC2
  map.idx <- match(rownames(estETOC2.m),rownames(data.m));
  rep.idx <- which(is.na(map.idx)==FALSE);
  print(paste("Number of represented epiTOC2 CpGs (max=163)=",length(rep.idx),sep=""))
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

age.v <- pheno$Age
epitoc_output <- epiTOC2(data.m = betas, age.v)
epitoc_tnsc<-as.data.frame(epitoc_output$tnsc)
colnames(epitoc_tnsc)<-"tnsc"
epitoc_irs<-as.data.frame(epitoc_output$irS)
colnames(epitoc_irs)<-"irs"
epitoc_pcgtage<-as.data.frame(epitoc_output$pcgtAge)
colnames(epitoc_pcgtage)<-"pcgtage"
