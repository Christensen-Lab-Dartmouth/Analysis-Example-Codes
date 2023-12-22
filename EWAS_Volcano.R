library(tidyr)
# this is the factor of interest
CaseStatus <- factor(Pheno$Case)
# this is the individual effect that we need to account for
Age <- Pheno$Age
Gender<-as.factor(Pheno$gender)

# use the above to create a design matrix
design <- model.matrix(~0+CaseStatus+Age+Gender)
colnames(design) <- c(levels(CaseStatus),"Age",levels(Gender)[-1])



library(limma)
fit1 <- lmFit(EWAS_Beta, design)

contMatrix1 <- makeContrasts(Tumor-Normal,levels = design)

fit2 <- contrasts.fit(fit1, contMatrix1)
fit2 <- eBayes(fit2)
delta_beta<-as.data.frame(fit2$coefficients)
colnames(delta_beta)[1]<-"Delta"
delta_beta$absDelta<-abs(delta_beta$Delta)
delta_beta<-rownames_to_column(delta_beta)
colnames(delta_beta)[1]<-"Name"

#delta_beta_filtered<-delta_beta %>% filter(absDelta >= 0.2)

# look at the numbers of DM CpGs at FDR < 0.01
summary(decideTests(fit2,adjust.method = "BH",p.value = 1e-10))

# get the table of results for the first contrast 
# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann450kSub <- ann450k[match(rownames(PhenoB),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]

DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
DMPs_delta<-inner_join(delta_beta,DMPs, by = "Name")

DMPsbuset<-DMPs_delta %>% filter(adj.P.Val < 1e-10 & absDelta >= 0.3)


write.xlsx(DMPsbuset,"CIMPHDMPsd0.3.xlsx")
write.xlsx(DMPs,"CIMPHDMPsComplete.xlsx")
table(DMPsbuset$logFC<0)#True=tumor-normal down=hypometh in tumor

DMPs_delta$log.P<-(-log10(DMPs_delta$adj.P.Val))
DMPs_delta$threshold<-ifelse(DMPs_delta$adj.P.Val<1e-10,"A","B")
DMPs_delta$threshold_D<-ifelse(DMPs_delta$adj.P.Val<1e-10 & DMPs_delta$absDelta >= 0.3,"A","B")


DMPs_delta %>% 
  group_by(threshold) %>% 
  summarise(count = sum(Delta >= 0.3), count2 = sum(Delta <= -0.3))

ggplot(DMPs_delta, aes(Delta, log.P)) +
  geom_point(aes(color = threshold_D),size = 1) +
  xlab("Delta Beta") + ylab("-log10(adjusted.P-value)")  +
  #scale_color_manual(values=c("blue","red"),name="No.CpGs",
  #labels=c("38467","16095"), guide = guide_legend(reverse=TRUE))+
  geom_vline(xintercept = 0.3, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = -0.3, color = "red", linetype = "dashed", size = 1) +
  geom_hline(yintercept = -log10(1e-10), color = "red", linetype = "dashed", size = 1)+
  scale_x_continuous(breaks = seq(-0.8,0.8,by = 0.2))+
  #scale_y_continuous(breaks = seq(0, 60, by= 10))+
  ggtitle("CIMP-H Tumor vs Normal")+
  theme(plot.title = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(colour="black", size = 20), legend.position = "none",
        #hjust = 1),
        axis.text.y=element_text(colour="black", size = 20),
        axis.title.x=element_text(colour="black", size = 20),
        axis.title.y=element_text(colour="black", size = 20),
        #legend.key = element_blank(),
        #legend.text = element_text (size = 17),
        #legend.title = element_text( size = 20)
  )