

library(pheatmap)
#Assume you want to color it by cell types
heat_annot <- data.frame(row.names = colnames(Beta),
                        CellType = Pheno$CellType)

pheatmap(
  Beta,
  annotation_names_col = F,
  show_rownames = TRUE, #CpGs 
  show_colnames = TRUE, #samples
  labels_col = Pheno$CellType,
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
