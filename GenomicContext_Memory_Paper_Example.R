bakcground_CpGs<-rownames(EWAS_results)
Bnv_Bmem<-EWAS_results %>% filter(P.value.adj.Bnv.minus.Bmem <0.05)
Bnv_Bmem$Coef.Bnv.minus.Bmem.abs<-abs(Bnv_Bmem$Coef.Bnv.minus.Bmem)
Bnv_Bmem<-Bnv_Bmem %>% filter(P.value.adj.Bnv.minus.Bmem <0.05 & Coef.Bnv.minus.Bmem.abs>0.2)
B_cpgs<-rownames(Bnv_Bmem)

CD4nv_CD4cm<-EWAS_results %>% filter(P.value.adj.CD4nv.minus.CD4CM <0.05)
CD4nv_CD4cm$Coef.CD4nv.minus.CD4CM.abs<-abs(CD4nv_CD4cm$Coef.CD4nv.minus.CD4CM)
CD4nv_CD4cm<-CD4nv_CD4cm %>% filter(P.value.adj.CD4nv.minus.CD4CM <0.05 & Coef.CD4nv.minus.CD4CM.abs>0.2)
CD4_cpgs<-rownames(CD4nv_CD4cm)

CD8nv_CD8cm<-EWAS_results %>% filter(P.value.adj.CD8nv.minus.CD8CM <0.05)
CD8nv_CD8cm$Coef.CD8nv.minus.CD8CM.abs<-abs(CD8nv_CD8cm$Coef.CD8nv.minus.CD8CM)
CD8nv_CD8cm<-CD8nv_CD8cm %>% filter(P.value.adj.CD8nv.minus.CD8CM <0.05 & Coef.CD8nv.minus.CD8CM.abs>0.2)
CD8_cpgs<-rownames(CD8nv_CD8cm)

library(reshape)
library(ggpubr)
library(epitools) #OR calculation
library(forestplot)
library(minfi)
annotation <- as.data.frame(getAnnotation("IlluminaHumanMethylationEPICanno.ilm10b2.hg19"))
Annotation <- annotation[annotation$Name %in% bakcground_CpGs,]
colnames(Annotation)[1]<-"chr"
Annotation$Relation_to_Island <- as.character(Annotation$Relation_to_Island )
# Add annotation collapsing shores and shelves for island context
Annotation$islandAdjust <- ifelse(Annotation$Relation_to_Island %in% c('N_Shelf','S_Shelf'),'Shelf',
                                  ifelse(Annotation$Relation_to_Island %in% c('N_Shore','S_Shore'),'Shore',Annotation$Relation_to_Island))

annotSub <- Annotation[,c('Name','chr','pos')]
epic.gr <- makeGRangesFromDataFrame(annotSub,keep.extra.columns = T,ignore.strand = T,seqnames.field = 'chr',start.field = 'pos',end.field = 'pos')
rm(annotSub)

library(genomation)
transcFeat <- readTranscriptFeatures('UCSC_hg19_refGene.bed.txt',up.flank = 2000, down.flank = 2000)
# The 'genomation' function annotates GRanges object as overlapping with promoter, exon, intron, or intergenic regions
epic.ann <- annotateWithGeneParts(epic.gr, transcFeat) 
membership <- epic.ann@members
# return %s for features 
getTargetAnnotationStats(epic.ann, percentage=TRUE, precedence=TRUE)

dim(membership)
membership2 <- membership
for(i in 1:nrow(membership2)){
  if(membership[i,1] == 1){membership2[i,] <- c(1,0,0)}
  if(membership[i,1] == 0 & membership[i,2] == 1){membership2[i,] <- c(0,1,0)}
  if(membership[i,1] == 0 & membership[i,2] == 0 & membership[i,3] == 1){membership2[1,] <- c(0,0,1)}
}

# add counts for transcription features to annotation 
Annotation$promoters <- membership2[,"prom"]
Annotation$exon <- membership2[,"exon"]
Annotation$intron <- membership2[,"intron"]
Annotation$intergenic <- 0
Annotation$intergenic[which(Annotation$promoters==0 & Annotation$exon==0 & Annotation$intron==0)] <- 1
# make sure TRUE totals to all probes 
table(Annotation$promoters)[[2]] + table(Annotation$exon)[[2]] + table(Annotation$intron)[[2]] + table(Annotation$intergenic)[[2]]

# Clean workspace
rm(transcFeat,epic.ann,membership,membership2)

Annotation$GC<-ifelse(Annotation$promoters == 1,"Promoters",ifelse(Annotation$exon == 1,"Exons",ifelse(Annotation$intron==1,"Introns","Intergenic")))

Annotation$P5EH<-ifelse(Annotation$Phantom5_Enhancers == "", "No", "Yes")

Annotation$BCPG<-ifelse(Annotation$Name %in% B_cpgs, 1, 0)
Annotation$CD4CPG<-ifelse(Annotation$Name %in% CD4_cpgs, 1, 0)
Annotation$CD8CPG<-ifelse(Annotation$Name %in% CD8_cpgs, 1, 0)
#B-------------------------------------------------------------------------------------------------
orTable2 <- as.data.frame(matrix(ncol = 5,nrow = 0))
colnames(orTable2) <- c('proportion','OR','lower','upper','pVal')

islandLevels <- c('Island','Shore','Shelf','OpenSea')
for (i in islandLevels) {
  prop <- table(Annotation$BCPG, Annotation$islandAdjust == i)['1','TRUE']/sum(Annotation$BCPG, na.rm = T)
  # Calculate odds ratio
  tempOR <- oddsratio.fisher(table(Annotation$BCPG,Annotation$islandAdjust == i))
  # Add values to output
  orTable2[i,] <- c(prop, tempOR$measure[2,],tempOR$p.value[2,'fisher.exact'])
}

orTable3 <- as.data.frame(matrix(ncol = 5,nrow = 0))
colnames(orTable3) <- c('proportion','OR','lower','upper','pVal')

islandLevels <- c('Promoters','Exons','Introns','Intergenic')
for (i in islandLevels) {
  prop <- table(Annotation$BCPG, Annotation$GC == i)['1','TRUE']/sum(Annotation$BCPG, na.rm = T)
  # Calculate odds ratio
  tempOR <- oddsratio.fisher(table(Annotation$BCPG,Annotation$GC == i))
  # Add values to output
  orTable3[i,] <- c(prop, tempOR$measure[2,],tempOR$p.value[2,'fisher.exact'])
}

orTable4 <- as.data.frame(matrix(ncol = 5,nrow = 0))
colnames(orTable4) <- c('proportion','OR','lower','upper','pVal')

islandLevels <- c('Yes','No')
for (i in islandLevels) {
  prop <- table(Annotation$BCPG, Annotation$P5EH == i)['1','TRUE']/sum(Annotation$BCPG, na.rm = T)
  # Calculate odds ratio
  tempOR <- oddsratio.fisher(table(Annotation$BCPG,Annotation$P5EH == i))
  # Add values to output
  orTable4[i,] <- c(prop, tempOR$measure[2,],tempOR$p.value[2,'fisher.exact'])
}
rownames(orTable4)[1]<-"Phantom5 Enhancers"

orTable<-rbind(orTable2,orTable3,orTable4)

orText <- cbind(c("Category","Relation to CpG Island","","","","Genomic Context","","","","Relation to Phantom5 Enhancers",""),
                c("Genomic Feature","CpG Island","Shore","Shelf","Open Sea","Promoters","Exons","Introns","Intergenic","Yes","No"),
                
                c('Odds Ratio (95% CI)',paste(format(round(orTable$OR,2),nsmall = 2),' (',format(round(orTable$lower,2),
                                                                                                 nsmall = 2),', ',
                                              format(round(orTable$upper,2),nsmall = 2),')',sep = '')))
# c('p-value',formatC(orTable$pVal,digits = 3)),
# c('Proportion of NR3C1 Probes (%)',paste0(format(round((orTable$proportion*100),2),nsmall = 2),"%"))


forestplot(labeltext = orText, 
           hrzl_lines = list("2" = gpar(lwd = 1, columns = 1:4, col = "#000044"),
                             "6" = gpar(lty = 2,columns = 1:4),
                             "10" = gpar(lty = 2,columns = 1:4),
                             "12" = gpar(lty = 1,columns = 1:4, col = "#000044")),
           #graph.pos = 2, 
           
           mean = c(NA,orTable$OR),
           lower = c(NA,orTable$lower),
           upper = c(NA,orTable$upper),
           #title="IPS DMCs NR3C1 Enrichment",
           #txt_gp=fpTxtGp(label=gpar(cex=1.3),ticks=gpar(cex=1.3)),
           is.summary=c(TRUE, FALSE, FALSE, FALSE,FALSE,FALSE,FALSE, FALSE, FALSE,FALSE,FALSE,FALSE,FALSE),
           txt_gp = fpTxtGp(summary = list(
             gpar(cex=1.2)),
             ticks=gpar(cex=1.2),xlab=gpar(cex=1.2)
           )
           , xlab = "Odds Ratio",
           
           # xticks = c(0,1,2,3,4,5,6,7,8,9),
           xlog = T,
           col=fpColors(box = "royalblue",
                        line = "darkblue")
           ,
           fn.ci_norm = fpDrawPointCI, pch = 16, cex=2, 
           zero = 1
           , 
           lineheight = "auto"
           ,
           boxsize=0.15
           , 
           colgap=unit(8,"mm")
)









#CD4-------------------------------------------------------------------------------------------------
orTable2 <- as.data.frame(matrix(ncol = 5,nrow = 0))
colnames(orTable2) <- c('proportion','OR','lower','upper','pVal')

islandLevels <- c('Island','Shore','Shelf','OpenSea')
for (i in islandLevels) {
  prop <- table(Annotation$CD4CPG, Annotation$islandAdjust == i)['1','TRUE']/sum(Annotation$CD4CPG, na.rm = T)
  # Calculate odds ratio
  tempOR <- oddsratio.fisher(table(Annotation$CD4CPG,Annotation$islandAdjust == i))
  # Add values to output
  orTable2[i,] <- c(prop, tempOR$measure[2,],tempOR$p.value[2,'fisher.exact'])
}

orTable3 <- as.data.frame(matrix(ncol = 5,nrow = 0))
colnames(orTable3) <- c('proportion','OR','lower','upper','pVal')

islandLevels <- c('Promoters','Exons','Introns','Intergenic')
for (i in islandLevels) {
  prop <- table(Annotation$CD4CPG, Annotation$GC == i)['1','TRUE']/sum(Annotation$CD4CPG, na.rm = T)
  # Calculate odds ratio
  tempOR <- oddsratio.fisher(table(Annotation$CD4CPG,Annotation$GC == i))
  # Add values to output
  orTable3[i,] <- c(prop, tempOR$measure[2,],tempOR$p.value[2,'fisher.exact'])
}

orTable4 <- as.data.frame(matrix(ncol = 5,nrow = 0))
colnames(orTable4) <- c('proportion','OR','lower','upper','pVal')

islandLevels <- c('Yes','No')
for (i in islandLevels) {
  prop <- table(Annotation$CD4CPG, Annotation$P5EH == i)['1','TRUE']/sum(Annotation$CD4CPG, na.rm = T)
  # Calculate odds ratio
  tempOR <- oddsratio.fisher(table(Annotation$CD4CPG,Annotation$P5EH == i))
  # Add values to output
  orTable4[i,] <- c(prop, tempOR$measure[2,],tempOR$p.value[2,'fisher.exact'])
}
rownames(orTable4)[1]<-"Phantom5 Enhancers"

orTable<-rbind(orTable2,orTable3,orTable4)

orText <- cbind(c("Category","Relation to CpG Island","","","","Genomic Context","","","","Relation to Phantom5 Enhancers",""),
                c("Genomic Feature","CpG Island","Shore","Shelf","Open Sea","Promoters","Exons","Introns","Intergenic","Yes","No"),
                
                c('Odds Ratio (95% CI)',paste(format(round(orTable$OR,2),nsmall = 2),' (',format(round(orTable$lower,2),
                                                                                                 nsmall = 2),', ',
                                              format(round(orTable$upper,2),nsmall = 2),')',sep = '')))
# c('p-value',formatC(orTable$pVal,digits = 3)),
# c('Proportion of NR3C1 Probes (%)',paste0(format(round((orTable$proportion*100),2),nsmall = 2),"%"))


forestplot(labeltext = orText, 
           hrzl_lines = list("2" = gpar(lwd = 1, columns = 1:4, col = "#000044"),
                             "6" = gpar(lty = 2,columns = 1:4),
                             "10" = gpar(lty = 2,columns = 1:4),
                             "12" = gpar(lty = 1,columns = 1:4, col = "#000044")),
           #graph.pos = 2, 
           
           mean = c(NA,orTable$OR),
           lower = c(NA,orTable$lower),
           upper = c(NA,orTable$upper),
           #title="IPS DMCs NR3C1 Enrichment",
           #txt_gp=fpTxtGp(label=gpar(cex=1.3),ticks=gpar(cex=1.3)),
           is.summary=c(TRUE, FALSE, FALSE, FALSE,FALSE,FALSE,FALSE, FALSE, FALSE,FALSE,FALSE,FALSE,FALSE),
           txt_gp = fpTxtGp(summary = list(
             gpar(cex=1.2)),
             ticks=gpar(cex=1.2),xlab=gpar(cex=1.2)
           )
           , xlab = "Odds Ratio",
           
           # xticks = c(0,1,2,3,4,5,6,7,8,9),
           xlog = T,
           col=fpColors(box = "royalblue",
                        line = "darkblue")
           ,
           fn.ci_norm = fpDrawPointCI, pch = 16, cex=2, 
           zero = 1
           , 
           lineheight = "auto"
           ,
           boxsize=0.15
           , 
           colgap=unit(8,"mm")
)









#CD8-------------------------------------------------------------------------------------------------
orTable2 <- as.data.frame(matrix(ncol = 5,nrow = 0))
colnames(orTable2) <- c('proportion','OR','lower','upper','pVal')

islandLevels <- c('Island','Shore','Shelf','OpenSea')
for (i in islandLevels) {
  prop <- table(Annotation$CD8CPG, Annotation$islandAdjust == i)['1','TRUE']/sum(Annotation$CD8CPG, na.rm = T)
  # Calculate odds ratio
  tempOR <- oddsratio.fisher(table(Annotation$CD8CPG,Annotation$islandAdjust == i))
  # Add values to output
  orTable2[i,] <- c(prop, tempOR$measure[2,],tempOR$p.value[2,'fisher.exact'])
}

orTable3 <- as.data.frame(matrix(ncol = 5,nrow = 0))
colnames(orTable3) <- c('proportion','OR','lower','upper','pVal')

islandLevels <- c('Promoters','Exons','Introns','Intergenic')
for (i in islandLevels) {
  prop <- table(Annotation$CD8CPG, Annotation$GC == i)['1','TRUE']/sum(Annotation$CD8CPG, na.rm = T)
  # Calculate odds ratio
  tempOR <- oddsratio.fisher(table(Annotation$CD8CPG,Annotation$GC == i))
  # Add values to output
  orTable3[i,] <- c(prop, tempOR$measure[2,],tempOR$p.value[2,'fisher.exact'])
}

orTable4 <- as.data.frame(matrix(ncol = 5,nrow = 0))
colnames(orTable4) <- c('proportion','OR','lower','upper','pVal')

islandLevels <- c('Yes','No')
for (i in islandLevels) {
  prop <- table(Annotation$CD8CPG, Annotation$P5EH == i)['1','TRUE']/sum(Annotation$CD8CPG, na.rm = T)
  # Calculate odds ratio
  tempOR <- oddsratio.fisher(table(Annotation$CD8CPG,Annotation$P5EH == i))
  # Add values to output
  orTable4[i,] <- c(prop, tempOR$measure[2,],tempOR$p.value[2,'fisher.exact'])
}
rownames(orTable4)[1]<-"Phantom5 Enhancers"

orTable<-rbind(orTable2,orTable3,orTable4)

orText <- cbind(c("Category","Relation to CpG Island","","","","Genomic Context","","","","Relation to Phantom5 Enhancers",""),
                c("Genomic Feature","CpG Island","Shore","Shelf","Open Sea","Promoters","Exons","Introns","Intergenic","Yes","No"),
                
                c('Odds Ratio (95% CI)',paste(format(round(orTable$OR,2),nsmall = 2),' (',format(round(orTable$lower,2),
                                                                                                 nsmall = 2),', ',
                                              format(round(orTable$upper,2),nsmall = 2),')',sep = '')))
# c('p-value',formatC(orTable$pVal,digits = 3)),
# c('Proportion of NR3C1 Probes (%)',paste0(format(round((orTable$proportion*100),2),nsmall = 2),"%"))


forestplot(labeltext = orText, 
           hrzl_lines = list("2" = gpar(lwd = 1, columns = 1:4, col = "#000044"),
                             "6" = gpar(lty = 2,columns = 1:4),
                             "10" = gpar(lty = 2,columns = 1:4),
                             "12" = gpar(lty = 1,columns = 1:4, col = "#000044")),
           #graph.pos = 2, 
           
           mean = c(NA,orTable$OR),
           lower = c(NA,orTable$lower),
           upper = c(NA,orTable$upper),
           #title="IPS DMCs NR3C1 Enrichment",
           #txt_gp=fpTxtGp(label=gpar(cex=1.3),ticks=gpar(cex=1.3)),
           is.summary=c(TRUE, FALSE, FALSE, FALSE,FALSE,FALSE,FALSE, FALSE, FALSE,FALSE,FALSE,FALSE,FALSE),
           txt_gp = fpTxtGp(summary = list(
             gpar(cex=1.2)),
             ticks=gpar(cex=1.2),xlab=gpar(cex=1.2)
           )
           , xlab = "Odds Ratio",
           
           # xticks = c(0,1,2,3,4,5,6,7,8,9),
           xlog = T,
           col=fpColors(box = "royalblue",
                        line = "darkblue")
           ,
           fn.ci_norm = fpDrawPointCI, pch = 16, cex=2, 
           zero = 1
           , 
           lineheight = "auto"
           ,
           boxsize=0.15
           , 
           colgap=unit(8,"mm")
)








