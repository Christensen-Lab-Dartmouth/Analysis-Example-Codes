CESC_Ango_Pred1<-CESC_Deconv3%>% dplyr::select("submitter_id", "tumor_stage", "vital_status", "days_to_last_follow_up", "days_to_death","Epithelial","Endothelial","Stromal")




CESC_Ango_Pred1$Survival_days<-ifelse(CESC_Ango_Pred1$days_to_last_follow_up == "NA",CESC_Ango_Pred1$days_to_death, 
                                      ifelse(CESC_Ango_Pred1$days_to_death == "NA",CESC_Ango_Pred1$days_to_last_follow_up,CESC_Ango_Pred1$days_to_death ))

CESC_Ango_Pred1<-CESC_Ango_Pred1 %>% filter(Survival_days >= 0 & Survival_days!="NA" & vital_status != "NA" & vital_status!="Not Reported")
CESC_Ango_Pred1$Survival_days<-as.numeric(CESC_Ango_Pred1$Survival_days)
CESC_Ango_Pred1$Angiogenic<-CESC_Ango_Pred1$Endothelial+CESC_Ango_Pred1$Epithelial+CESC_Ango_Pred1$Stromal

ggplot(CESC_Ango_Pred1, aes(x = vital_status, y = Angiogenic))+ 
  xlab("") +
  #ylab("CSC (%)") + #ylim(35,82)+
  geom_boxplot(outlier.size=0.5, fatten = 1)+
  stat_compare_means(method = "wilcox.test", label = "p", size = 5, comparisons = list(c("Alive","Dead") ))+
  
  theme( axis.text = element_text( size = 10 ),
         axis.text.x = element_text(size = 15, face = "bold" ),
         axis.title = element_text( size = 15 , face="bold"))  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #legend.position = "none"
  )
#Tumor
library(pROC)
library(survival) 
library(survminer)
model <- survfit(Surv(Survival_days,vital_status == "Dead") ~ Angiogenic>median(CESC_Ango_Pred1$Angiogenic,na.rm = T),
                 data = CESC_Ango_Pred1)
summary(model)
ggsurvplot(model,  
           pval = TRUE, #conf.int = TRUE,
           risk.table = F,ggtheme = theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")),
           legend.labs =  c("Angiogenic% \u2264 3.91%", "Angiogenic% > 3.91%"), 
           title = "CESC Angiogenic% and survival")
#---------------------------------------------------------------------------------------
library(gtools)
CESC_Ango_Pred1$Endo_quant<-quantcut(CESC_Ango_Pred1$Angiogenic, q = 4)
model <- survfit(Surv(Survival_days,vital_status == "Dead") ~ Endo_quant,
                 data = CESC_Ango_Pred1)
summary(model)
ggsurvplot(model,  
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE,ggtheme = theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")),
           # legend.labs =  c("Tumor% \u2264 78.1%", "Tumor% > 78.1%"), 
           title = "CESC Tumor% and survival")
#--------------------------------------------------------------------------------------------
ggplot(CESC_Ango_Pred1, aes(x = vital_status, y = Endothelial))+ 
  xlab("") +
  #ylab("CSC (%)") + #ylim(35,82)+
  geom_boxplot(outlier.size=0.5, fatten = 1)+
  stat_compare_means(method = "wilcox.test", label = "p", size = 5, comparisons = list(c("Alive","Dead") ))+
  
  theme( axis.text = element_text( size = 10 ),
         axis.text.x = element_text(size = 15, face = "bold" ),
         axis.title = element_text( size = 15 , face="bold"))  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #legend.position = "none"
  )
#Tumor
library(pROC)
library(survival) 
library(survminer)
model <- survfit(Surv(Survival_days,vital_status == "Dead") ~ Endothelial>median(CESC_Ango_Pred1$Endothelial,na.rm = T),
                 data = CESC_Ango_Pred1)
summary(model)
ggsurvplot(model,  
           pval = TRUE, #conf.int = TRUE,
           risk.table = F,ggtheme = theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black")),
           legend.labs =  c("Endothelial% \u2264 3.91%", "Endothelial% > 3.91%"), 
           title = "CESC Endothelial% and survival")
#---------------------------------------------------------------------------------------
library(gtools)
CESC_Ango_Pred1$Endo_quant<-quantcut(CESC_Ango_Pred1$Endothelial, q = 3)
model <- survfit(Surv(Survival_days,vital_status == "Dead") ~ Endo_quant,
                 data = CESC_Ango_Pred1)
summary(model)
ggsurvplot(model,  
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE,ggtheme = theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")),
           # legend.labs =  c("Tumor% \u2264 78.1%", "Tumor% > 78.1%"), 
           title = "CESC Tumor% and survival")
#--------------------------------------------------------------------------------------------
ggplot(CESC_Ango_Pred1, aes(x = vital_status, y = Epithelial))+ 
  xlab("") +
  #ylab("CSC (%)") + #ylim(35,82)+
  geom_boxplot(outlier.size=0.5, fatten = 1)+
  stat_compare_means(method = "wilcox.test", label = "p", size = 5, comparisons = list(c("Alive","Dead") ))+
  
  theme( axis.text = element_text( size = 10 ),
         axis.text.x = element_text(size = 15, face = "bold" ),
         axis.title = element_text( size = 15 , face="bold"))  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #legend.position = "none"
  )
#Tumor
library(pROC)
library(survival) 
library(survminer)
model <- survfit(Surv(Survival_days,vital_status == "Dead") ~ Epithelial>median(CESC_Ango_Pred1$Epithelial,na.rm = T),
                 data = CESC_Ango_Pred1)
summary(model)
ggsurvplot(model,  
           pval = TRUE, #conf.int = TRUE,
           risk.table = F,ggtheme = theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black")),
           legend.labs =  c("Epithelial% \u2264 3.91%", "Epithelial% > 3.91%"), 
           title = "CESC Epithelial% and survival")
#---------------------------------------------------------------------------------------
library(gtools)
CESC_Ango_Pred1$Epithelial_Tritile<-quantcut(CESC_Ango_Pred1$Epithelial, q = 3)
model <- survfit(Surv(Survival_days,vital_status == "Dead") ~ Epithelial_Tritile,
                 data = CESC_Ango_Pred1)
summary(model)
ggsurvplot(model,  
           pval = TRUE, #conf.int = TRUE,
           risk.table = FALSE,ggtheme = theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")),
           # legend.labs =  c("Tumor% \u2264 78.1%", "Tumor% > 78.1%"), 
           title = "CESC Epithelial% and survival")
#--------------------------------------------------------------------------------------------
ggplot(CESC_Ango_Pred1, aes(x = vital_status, y = Stromal))+ 
  xlab("") +
  #ylab("CSC (%)") + #ylim(35,82)+
  geom_boxplot(outlier.size=0.5, fatten = 1)+
  stat_compare_means(method = "wilcox.test", label = "p", size = 5, comparisons = list(c("Alive","Dead") ))+
  
  theme( axis.text = element_text( size = 10 ),
         axis.text.x = element_text(size = 15, face = "bold" ),
         axis.title = element_text( size = 15 , face="bold"))  + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        #legend.position = "none"
  )
#Tumor
library(pROC)
library(survival) 
library(survminer)
model <- survfit(Surv(Survival_days,vital_status == "Dead") ~ Stromal>median(CESC_Ango_Pred1$Stromal,na.rm = T),
                 data = CESC_Ango_Pred1)
summary(model)
ggsurvplot(model,  
           pval = TRUE, #conf.int = TRUE,
           risk.table = F,ggtheme = theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black")),
           legend.labs =  c("Stromal% \u2264 4.17%", "Stromal% > 4.17%"), 
           title = "CESC Stromal% and survival")
#---------------------------------------------------------------------------------------
library(gtools)
CESC_Ango_Pred1$Endo_quant<-quantcut(CESC_Ango_Pred1$Stromal, q = 3)
model <- survfit(Surv(Survival_days,vital_status == "Dead") ~ Endo_quant,
                 data = CESC_Ango_Pred1)
summary(model)
ggsurvplot(model,  
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE,ggtheme = theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black")),
           # legend.labs =  c("Tumor% \u2264 78.1%", "Tumor% > 78.1%"), 
           title = "CESC Tumor% and survival")
