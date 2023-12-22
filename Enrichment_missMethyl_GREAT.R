
library(missMethyl)
library(methylGSA)
library(tidyr)
library(limma)
library(ddpcr)
library(dplyr)
library(tibble)
library(openxlsx)

#missMethy
asGO<-gometh(
  AsCpGs, #interested CpGs
  all.cpg = allCpGs, #Background CpGs
  collection = "GO",
  array.type = "450K",
  plot.bias = FALSE,
  prior.prob = TRUE,
  anno = NULL,
  equiv.cpg = TRUE,
  fract.counts = TRUE
)

asGO1 <- asGO %>% filter(N > 5 & N < 2000)
asGO1 <- asGO1[order(asGO1$FDR),] %>% filter(FDR< 0.05)


#----------------------------------------------------------------------------------------------------------------------

#GREAT

library(rtracklayer)
gran<-as.data.frame(hm450.manifest)
subgran<-gran[rownames(gran) %in% IPS_all_CpGs,]
subgran <- rownames_to_column(subgran)
colnames(subgran)[1]<-"name"
subgranges<- subgran %>% dplyr::select(seqnames, start, end, name) 
subgranges1<-GRanges(subgranges)
export(subgranges1, "background.bed", format = "bed")

subgran<-gran[rownames(gran) %in% AsCpGs,]
subgran <- rownames_to_column(subgran)
colnames(subgran)[1]<-"name"
subgranges<- subgran %>% dplyr::select(seqnames, start, end, name) 
subgranges1<-GRanges(subgranges)
export(subgranges1, "Astrocytes.bed", format = "bed")



# get tables for each GO set ordered by variable of interest 
get_GO_tables <- function(results, GO_table_to_get, order_by){
  tb = getEnrichmentTables(results)
  # restrict to which GP processes are desired 
  tb <- tb[[GO_table_to_get]]
  # calculate FDR from raw P-values 
  tb$Hyper_FDR_Q_val <- p.adjust(tb$Hyper_Raw_PValue, method = "BH")
  # restrict to GO terms w/ FDR<0.01
  tb <- tb[tb$Hyper_FDR_Q_val<0.01,]
  # order results by FDR or enrichment value 
  if(order_by == "qval") 
    tb <- tb[order(tb$Hyper_FDR_Q_val, decreasing = FALSE),]
  if(order_by == "enrichment") 
    tb <- tb[order(tb$Hyper_Fold_Enrichment, decreasing = TRUE),]
  # return processed table 
  tb
}


library(rGREAT)
library(genomation)
background <- readBed('background.bed') 

Ast <- readBed('Astrocytes.bed') 
# submit GREAT job (specifying background set as 450K probes that passed QC)
job = submitGreatJob(gr = Ast, bg = background, species = "hg19",request_interval = 0)

# GO Biological Process
tb_bp_qval_ast <- get_GO_tables(job, "GO Biological Process", "qval")

# GO Molecular Function 
tb_mf_qval_ast <- get_GO_tables(job, "GO Molecular Function", "qval")
# GO Cellular Component
tb_cc_qval_ast <- get_GO_tables(job, "GO Cellular Component", "qval")
#combined
ast_GREAT<-rbind(tb_bp_qval_ast,tb_mf_qval_ast,tb_cc_qval_ast)
ast_GREAT <- ast_GREAT[order(ast_GREAT$Hyper_FDR_Q_val), ] %>% filter(Total_Genes_Annotated > 5 & Total_Genes_Annotated < 2000)
ast_GREAT$Category<-ifelse(ast_GREAT$ID %in% tb_bp_qval_ast$ID, "BP",ifelse(ast_GREAT$ID %in% tb_mf_qval_ast$ID,"MF","CC"))
ast_GREAT_search<-dplyr::filter(ast_GREAT , grepl("ast",name))

