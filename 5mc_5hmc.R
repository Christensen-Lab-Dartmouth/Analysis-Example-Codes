workdir<-"."
IDATprefixes <- searchIDATprefixes(workdir, recursive=TRUE)
register(FutureParam())
cl <- parallel::makeCluster(detectCores(), type = "SOCK")
plan(cluster, workers = cl,.cleanup = TRUE)
p<-FutureParam()
tmp2 <- openSesame(IDATprefixes,quality.mask = TRUE,
                   nondetection.mask = TRUE,
                   correct.switch = TRUE,
                   mask.use.tcga = FALSE,
                   pval.threshold = 0.05,
                   pval.method = "pOOBAH",
                   sum.TypeI = TRUE, 
                   platform = "EPIC", 
                   BPPARAM = p,
                   what="sigset")
parallel::stopCluster(cl)
plan(sequential)
gc()
closeAllConnections()

register(FutureParam())
cl <- parallel::makeCluster(detectCores(), type = "SOCK")
plan(cluster, workers = cl,.cleanup = TRUE)
p<-FutureParam()
tmp3<-SigSetsToRGChannelSet(tmp2, BPPARAM = p)
parallel::stopCluster(cl)
plan(sequential)
gc()
closeAllConnections()

BS_RGCS<-tmp3


pheno<-read.csv("1_SampleAnnotation.csv")
plotCtrl2(BS_RGCS)
densityPlot(BS_RGCS, xlab = "BS Beta")
#---------------------------------------------------------------------------------------------------------------
workdir<-"."
IDATprefixes <- searchIDATprefixes(workdir, recursive=TRUE)
register(FutureParam())
cl <- parallel::makeCluster(detectCores(), type = "SOCK")
plan(cluster, workers = cl,.cleanup = TRUE)
p<-FutureParam()
tmp2 <- openSesame(IDATprefixes,quality.mask = TRUE,
                   nondetection.mask = TRUE,
                   correct.switch = TRUE,
                   mask.use.tcga = FALSE,
                   pval.threshold = 0.05,
                   pval.method = "pOOBAH",
                   sum.TypeI = TRUE, 
                   platform = "EPIC", 
                   BPPARAM = p,
                   what="sigset")
parallel::stopCluster(cl)
plan(sequential)
gc()
closeAllConnections()

register(FutureParam())
cl <- parallel::makeCluster(detectCores(), type = "SOCK")
plan(cluster, workers = cl,.cleanup = TRUE)
p<-FutureParam()
tmp3<-SigSetsToRGChannelSet(tmp2, BPPARAM = p)
parallel::stopCluster(cl)
plan(sequential)
gc()
closeAllConnections()

oxBS_RGCS<-tmp3

plotCtrl(oxBS_RGCS)
densityPlot(oxBS_RGCS, xlab = "oxBS Beta")


#-------------------------------------------------------------------------------------------------------------
All_beta<-openSesame(".", platform = "EPIC")
sum(is.na(All_beta))
densityPlot(All_beta)
pheno<-read.csv("1_SampleAnnotation.csv")


pheno1<-pheno[order(pheno$ï..ChipID),]
identical(pheno1$ï..ChipID,colnames(All_beta))

densityPlot(All_beta, main="After QC", sampGroups = pheno1$Conversion,
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(pheno1$Conversion)) ,text.col=brewer.pal(8,"Dark2"))




All_beta_oxBS<-All_beta[,c(1,3,4,7,10,11,13,15)]
All_beta_BS<-All_beta[,c(2,5,6,8,9,12,14,16)]


colnames(All_beta_oxBS)<-pheno1[c(1,3,4,7,10,11,13,15),2]
colnames(All_beta_BS)<-pheno1[c(2,5,6,8,9,12,14,16),2]


All_beta_oxBS1<-All_beta_oxBS[,order(colnames(All_beta_oxBS))]
All_beta_BS1<-All_beta_BS[,order(colnames(All_beta_BS))]

identical(colnames(All_beta_BS1),colnames(All_beta_oxBS1))
identical(rownames(All_beta_BS1),rownames(All_beta_oxBS1))


NA_ID<-All_beta_BS1+All_beta_oxBS1



workdir<-"."
IDATprefixes <- searchIDATprefixes(workdir, recursive=TRUE)
register(FutureParam())
cl <- parallel::makeCluster(detectCores(), type = "SOCK")
plan(cluster, workers = cl,.cleanup = TRUE)
p<-FutureParam()
tmp2 <- openSesame(IDATprefixes,quality.mask = TRUE,
                   nondetection.mask = TRUE,
                   correct.switch = TRUE,
                   mask.use.tcga = FALSE,
                   pval.threshold = 0.05,
                   pval.method = "pOOBAH",
                   sum.TypeI = TRUE, 
                   platform = "EPIC", 
                   BPPARAM = p,
                   what="sigset")
parallel::stopCluster(cl)
plan(sequential)
gc()
closeAllConnections()

register(FutureParam())
cl <- parallel::makeCluster(detectCores(), type = "SOCK")
plan(cluster, workers = cl,.cleanup = TRUE)
p<-FutureParam()
tmp3<-SigSetsToRGChannelSet(tmp2, BPPARAM = p)
parallel::stopCluster(cl)
plan(sequential)
gc()
closeAllConnections()

RGCS<-tmp3

plotCtrl(RGCS)




MSet <- preprocessRaw(RGCS) 
RGCS_beta<-getBeta(MSet)


densityPlot(RGCS_beta, main="Before QC", sampGroups = pheno1$Conversion,
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(pheno1$Conversion)) ,text.col=brewer.pal(8,"Dark2"))




Meth<-getMeth(MSet)
Unmeth<-getUnmeth(MSet)
Total_Signal<-Meth+Unmeth

identical(colnames(RGCS_beta),pheno1$ï..ChipID)
identical(colnames(Total_Signal),pheno1$ï..ChipID)

RGCS_beta_oxBS<-RGCS_beta[,c(1,3,4,7,10,11,13,15)]
RGCS_beta_BS<-RGCS_beta[,c(2,5,6,8,9,12,14,16)]

colnames(RGCS_beta_oxBS)<-pheno1[c(1,3,4,7,10,11,13,15),2]
colnames(RGCS_beta_BS)<-pheno1[c(2,5,6,8,9,12,14,16),2]

Signal_oxBS<-Total_Signal[,c(1,3,4,7,10,11,13,15)]
Signal_BS<-Total_Signal[,c(2,5,6,8,9,12,14,16)]

colnames(Signal_oxBS)<-pheno1[c(1,3,4,7,10,11,13,15),2]
colnames(Signal_BS)<-pheno1[c(2,5,6,8,9,12,14,16),2]


RGCS_beta_oxBS1<-RGCS_beta_oxBS[,order(colnames(RGCS_beta_oxBS))]
RGCS_beta_BS1<-RGCS_beta_BS[,order(colnames(RGCS_beta_BS))]

Signal_oxBS1<-Signal_oxBS[,order(colnames(Signal_oxBS))]
Signal_BS1<-Signal_BS[,order(colnames(Signal_BS))]


hmc<-oxBS.MLE(RGCS_beta_BS1,RGCS_beta_oxBS1,Signal_BS1,Signal_oxBS1)
densityPlot(hmc$`5mC`,main = "5mC before QC")
densityPlot(hmc$`5hmC`, main = "5hmC before QC")

mc5<-hmc$`5mC`
mhc5<-hmc$`5hmC`

mc5<-mc5[order(rownames(mc5)),]
mhc5<-mhc5[order(rownames(mhc5)),]
NA_ID<-NA_ID[order(rownames(NA_ID)),]
identical(colnames(NA_ID),colnames(mc5))
identical(rownames(NA_ID),rownames(mhc5))

rowID<-intersect(rownames(NA_ID),rownames(mc5))

mc5<-mc5[rownames(mc5) %in% rowID,]
mhc5<-mhc5[rownames(mhc5) %in% rowID,]
NA_ID<-NA_ID[rownames(NA_ID) %in% rowID,]

NA_coord<-as.data.frame(which(is.na(NA_ID), arr.ind=TRUE))
for (i in 1:nrow(NA_coord)) {
  
    x<-NA_coord[i,1]
    y<-NA_coord[i,2]
    mc5[x,y]<-NA
  
  
}

for (i in 1:nrow(NA_coord)) {
  
  x<-NA_coord[i,1]
  y<-NA_coord[i,2]
  mhc5[x,y]<-NA
  
  
}

densityPlot(mc5,main = "5mC After QC")
densityPlot(mhc5, main = "5hmC After QC")

unfiltered_5mC<-hmc$`5mC`
unfiltered_5mC<-unfiltered_5mC[order(rownames(unfiltered_5mC)),]
unfiltered_5hmC<-hmc$`5hmC`
unfiltered_5hmC<-unfiltered_5hmC[order(rownames(unfiltered_5hmC)),]
filtered_5mC<-mc5
filtered_5hmC<-mhc5

save(unfiltered_5mC,unfiltered_5hmC,filtered_5mC,filtered_5hmC, file = "Pilot_Kidney_5mc_5hmC.RDATA")
#--------------------------------------------------------------------------------------------------------------
boxplot(CCRCC_pheno$GCT~CCRCC_pheno$OxBS)
boxplot(CCRCC_pheno$frac_na~CCRCC_pheno$OxBS)
boxplot(CCRCC_pheno$mean_intensity~CCRCC_pheno$OxBS)
boxplot(CCRCC_pheno$mean_beta~CCRCC_pheno$OxBS)

boxplot(CCRCC_pOOBHA)
identical(colnames(CCRCC_pOOBHA),CCRCC_pheno$Sample_ID)
CCRCC_pOOBHA1<-CCRCC_pOOBHA
colnames(CCRCC_pOOBHA1)<-CCRCC_pheno$OxBS
boxplot(CCRCC_pOOBHA1)

pval<-detectionP(CCRCC_rgset)
pval[1:10,1:10]
CCRCC_pOOBHA[1:10,1:10]
boxplot(pval)

pval["cg12950382",]
CCRCC_pOOBHA["cg12950382",]

identical(colnames(pval),colnames(CCRCC_pOOBHA))
identical(colnames(pval),CCRCC_pheno$Sample_ID)
oxBS_ID<-CCRCC_pheno[CCRCC_pheno$OxBS == "oxBS1", "Sample_ID"]
BS_ID<-CCRCC_pheno[CCRCC_pheno$OxBS == "BS1", "Sample_ID"]

BS_pval<-CCRCC_pOOBHA[,BS_ID]
oxBS_pval<-pval[,oxBS_ID]
BS_pval<-BS_pval[rownames(BS_pval) %in% rownames(oxBS_pval),]
combined_pval<-cbind(BS_pval,oxBS_pval)
combined_pval<-combined_pval[,order(colnames(combined_pval))]
CCRCC_pheno<-CCRCC_pheno[order(CCRCC_pheno$Sample_ID),]
identical(colnames(combined_pval),CCRCC_pheno$Sample_ID)
combined_pval1<-combined_pval
colnames(combined_pval1)<-CCRCC_pheno$OxBS
boxplot(combined_pval1)

MSet <- preprocessRaw(CCRCC_rgset) 
RGCS_beta<-getBeta(MSet)
identical(colnames(RGCS_beta),colnames(combined_pval))

RGCS_beta<-RGCS_beta[order(rownames(RGCS_beta)),]
combined_pval<-combined_pval[order(rownames(combined_pval)),]
identical(rownames(RGCS_beta),rownames(combined_pval))

for (i in 1:nrow(combined_pval)) {
  for (j in 1:ncol(combined_pval)) {
    
    if (combined_pval[i,j]<0.05) {
      RGCS_beta[i,j] = RGCS_beta[i,j] 
    } else {
      RGCS_beta[i,j] = NA
    }
    
  }
}


na_frac<-as.data.frame(colSums(is.na(RGCS_beta))/nrow(RGCS_beta))
colnames(na_frac)<-"NA Fraction"
identical(rownames(na_frac),CCRCC_pheno$Sample_ID)
na_frac$oxBS<-CCRCC_pheno$OxBS
boxplot(na_frac$`NA Fraction`~na_frac$oxBS)


#-----------------------------------------------------------------------------------------
MSet <- preprocessRaw(CCRCC_rgset) 
RGCS_beta<-getBeta(MSet)
CCRCC_pOOBHA2<-CCRCC_pOOBHA[rownames(CCRCC_pOOBHA) %in% rownames(RGCS_beta),]
identical(colnames(RGCS_beta),colnames(CCRCC_pOOBHA2))

RGCS_beta<-RGCS_beta[order(rownames(RGCS_beta)),]
CCRCC_pOOBHA2<-CCRCC_pOOBHA2[order(rownames(CCRCC_pOOBHA2)),]
identical(rownames(RGCS_beta),rownames(CCRCC_pOOBHA2))

for (i in 1:nrow(CCRCC_pOOBHA2)) {
  for (j in 1:ncol(CCRCC_pOOBHA2)) {
    
    if (CCRCC_pOOBHA2[i,j]<0.05) {
      RGCS_beta[i,j] = RGCS_beta[i,j] 
    } else {
      RGCS_beta[i,j] = NA
    }
    
  }
}


na_frac2<-as.data.frame(colSums(is.na(RGCS_beta))/nrow(RGCS_beta))
colnames(na_frac2)<-"NA Fraction"
identical(rownames(na_frac2),CCRCC_pheno$Sample_ID)
na_frac2$oxBS<-CCRCC_pheno$OxBS
boxplot(na_frac2$`NA Fraction`~na_frac2$oxBS)


#-----------------------------------------------------------------------------------------------------
oxBS_ID<-CCRCC_pheno[CCRCC_pheno$OxBS == "oxBS1", "Sample_ID"]
BS_ID<-CCRCC_pheno[CCRCC_pheno$OxBS == "BS1", "Sample_ID"]

BS_pval<-CCRCC_pOOBHA[,BS_ID]
oxBS_pval<-CCRCC_pOOBHA[,oxBS_ID]

identical(rownames(BS_pval),rownames(oxBS_pval))
BS_pval_sig_pc<-matrix(NA,nrow(BS_pval),1)

for (i in 1:nrow(BS_pval)) {
BS_pval_sig_pc[i,1]<-sum(BS_pval[i,]<0.05)/ncol(BS_pval)
}
rownames(BS_pval_sig_pc)<-rownames(BS_pval)
colnames(BS_pval_sig_pc)<-"BS Significant proportion"


oxBS_pval_sig_pc<-matrix(NA,nrow(oxBS_pval),1)

for (i in 1:nrow(oxBS_pval)) {
  oxBS_pval_sig_pc[i,1]<-sum(oxBS_pval[i,]<0.05)/ncol(oxBS_pval)
}
rownames(oxBS_pval_sig_pc)<-rownames(oxBS_pval)
colnames(oxBS_pval_sig_pc)<-"oxBS Significant proportion"
pval_sig_pc<-as.data.frame(cbind(BS_pval_sig_pc,oxBS_pval_sig_pc))
pval_sig_pc$diff<-pval_sig_pc[,1]-pval_sig_pc[,2]

rm(list = ls()[!ls() %in% c("pval_sig_pc")])
#---------------------------------------------------------------------------------------------
CCRCC_pheno<-CCRCC_pheno[order(CCRCC_pheno$Sample_Name),]
CCRCC_pheno[1,2]<-"RRC165_oxBS1"
CCRCC_pOOBHA<-CCRCC_pOOBHA[,CCRCC_pheno$Sample_ID]
identical(CCRCC_pheno$Sample_ID, colnames(CCRCC_pOOBHA))
boxplot(CCRCC_pOOBHA)




CCRCC_pOOBHA1<-CCRCC_pOOBHA
colnames(CCRCC_pOOBHA1)<-CCRCC_pheno$OxBS
boxplot(CCRCC_pOOBHA1)

CCRCC_pheno1<-CCRCC_pheno %>% dplyr::select(Sample_Name,Subject,OxBS)

paired.pval<-matrix(NA,nrow(CCRCC_pOOBHA),1)

for (i in 1:nrow(CCRCC_pOOBHA)) {
  r1<-CCRCC_pOOBHA[i,]
  CCRCC_pheno2<-cbind(r1,CCRCC_pheno1)
  CCRCC_pheno2<-CCRCC_pheno2[order(CCRCC_pheno2$OxBS),]
  x<-wilcox.test(r1 ~ OxBS, data = CCRCC_pheno2, paired = TRUE, alternative = "less")$p.value
  paired.pval[i,1]<-x
}
hist(CCRCC_pOOBHA)
rownames(paired.pval)<-rownames(CCRCC_pOOBHA)
paired.pval<-as.data.frame(paired.pval)
paired.pval$FDR<-p.adjust(paired.pval$V1, method = "fdr", n = nrow(paired.pval))
sum(paired.pval$V1<0.05)
CCRCC_paired_pval<-paired.pval

sum(CCRCC_paired_pval_PooBHA$FDR<0.05)
#--------------------------------------------------------------------------------
CCRCC_pheno<-CCRCC_pheno[order(CCRCC_pheno$OxBS),]
CCRCC_pheno<-CCRCC_pheno[order(CCRCC_pheno$Subject),]

bs_pheno<-CCRCC_pheno  %>% select(Subject,OxBS,frac_na) %>% filter(OxBS == "BS1")
oxBS_pheno<-CCRCC_pheno %>% select(Subject,OxBS,frac_na)  %>% filter(OxBS == "oxBS1")
identical(bs_pheno$Subject,oxBS_pheno$Subject)
paired_pheno<-rbind(bs_pheno,oxBS_pheno)
ggpaired(paired_pheno, x = "OxBS", y = "frac_na",
         color = "OxBS", line.color = "gray", line.size = 0.4,
         palette = "jco")+
  stat_compare_means(paired = TRUE)+xlab("")+ylab("Insignificant probe (NA) %") + ggtitle("CCRCC")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))

oxBS_ID<-rownames(CCRCC_pheno[CCRCC_pheno$OxBS== "oxBS1",])
BS_ID<-rownames(CCRCC_pheno[CCRCC_pheno$OxBS == "BS1",])

BS_pval<-CCRCC_pOOBHA[,BS_ID]
oxBS_pval<-CCRCC_pOOBHA[,oxBS_ID]

BS_pval_stack<-as.data.frame(stack(BS_pval))
BS_pval_stack<-BS_pval_stack[,2:3]
BS_pval_stack[,1]<-"BS"

oxBS_pval_stack<-as.data.frame(stack(oxBS_pval))
oxBS_pval_stack<-oxBS_pval_stack[,2:3]
oxBS_pval_stack[,1]<-"oxBS"

pval_stack<-rbind(BS_pval_stack,oxBS_pval_stack)
#pval_stack$col<- factor(pval_stack$col , levels=c("BS","oxBS"))
ggplot(data=pval_stack, aes(x=value, group=col, fill=col)) +
  geom_density(adjust=1.5) +scale_x_continuous(breaks=seq(0, 0.1, 0.02), limits=c(0, 0.1))+xlab("pOOBHA")+
  ggtitle("ccRCC")+
  #theme_ipsum() +
  facet_wrap(~col) +
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))

