EWAS_results<-read.csv("EWAS_Memory.csv", row.names = 1)

CD4cm_CD4em<-EWAS_results %>% filter(P.value.adj.CD4EM.minus.CD4CM <0.05)
CD4cm_CD4em$Coef.CD4EM.minus.CD4CM.abs<-abs(CD4cm_CD4em$Coef.CD4EM.minus.CD4CM)
CD4cm_CD4em<-CD4cm_CD4em[order(CD4cm_CD4em$Coef.CD4EM.minus.CD4CM.abs, decreasing = T),]
CD4cm_CD4em_20<-rownames(CD4cm_CD4em[1:20,])

CD8cm_CD8em<-EWAS_results %>% filter(P.value.adj.CD8EM.minus.CD8CM <0.05)
CD8cm_CD8em$Coef.CD8EM.minus.CD8CM.abs<-abs(CD8cm_CD8em$Coef.CD8EM.minus.CD8CM)
CD8cm_CD8em<-CD8cm_CD8em[order(CD8cm_CD8em$Coef.CD8EM.minus.CD8CM.abs, decreasing = T),]
CD8cm_CD8em_20<-rownames(CD8cm_CD8em[1:20,])

#heatmaps------------------------------------------------------------------------
CD4_pheno<-pheno %>% filter(CellType %in% c("CD4EM","CD4CM"))
CD4_pheno$CellType <- factor(CD4_pheno$CellType, levels = c("CD4CM","CD4EM"))

#CD4_pheno<-CD4_pheno[order(CD4_pheno$CellType),]

CD4beta<-as.matrix(beta[CD4cm_CD4em_20,rownames(CD4_pheno)])



identical(colnames(CD4beta),rownames(CD4_pheno))
heat_annot <- data.frame(row.names = colnames(CD4beta),
                         CellType = CD4_pheno$CellType)


library(pheatmap)
pheatmap(
  CD4beta,
  annotation_names_col = F,
  show_rownames = TRUE, #CpGs 
  show_colnames = TRUE, #samples
  labels_col = CD4_pheno$CellType,
  annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  #clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  #clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

CD4cmem_annot<-annot[CD4cm_CD4em_20,]

#------------------------------------------------------------------------------
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
CD8_pheno<-pheno %>% filter(CellType %in% c("CD8EM","CD8CM"))
CD8_pheno$CellType <- factor(CD8_pheno$CellType, levels = c("CD8CM","CD8EM"))

#CD8_pheno<-CD8_pheno[order(CD8_pheno$CellType),]

CD8beta<-as.matrix(beta[CD8cm_CD8em_20,rownames(CD8_pheno)])



identical(colnames(CD8beta),rownames(CD8_pheno))
heat_annot <- data.frame(row.names = colnames(CD8beta),
                         CellType = CD8_pheno$CellType)


library(pheatmap)
pheatmap(
  CD8beta,
  annotation_names_col = F,
  show_rownames = TRUE, #CpGs 
  show_colnames = TRUE, #samples
  labels_col = CD8_pheno$CellType,
  annotation_col = heat_annot,
  #annotation_colors = ann_colors,
  color = colorRampPalette(c("yellow","black","blue"))(128),
  clustering_distance_rows = "manhattan",
  #clustering_distance_cols = "manhattan",
  #cluster_rows=FALSE,
  clustering_method = "average",
  clustering_callback = callback,
  border_color = NA,
  fontsize = 10
)

CD8cmem_annot<-annot[CD8cm_CD8em_20,]




