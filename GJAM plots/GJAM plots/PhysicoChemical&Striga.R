library(gplots)
library(gjam)
library(tidyr)
library(ggplot2)
library(dplyr)
setwd("C:/Users/small/Documents/Brady lab/Soil Paper/Marcio")
load("C:/Users/small/Documents/Brady lab/Soil Paper/Marcio/GJAMResults_MLNew.RData")

SigVarsSeed<- as.character(SigVars.seed$TaxVar) #attempt at making SigVarsSeed as one wasn't provided


#Looking for physical-chemical factor that correlates with SEEDBANK--------------------
(SigVarsSeed.Soil <- SigVarsSeed[!grepl("ASV|Stri", SigVarsSeed)]) #grab what does not have ASV or Striga in name

distance= dist(t(netTab.seed[,SigVarsSeed.Soil]), method ="euclidian")    
hcluster = hclust(distance, method ="complete")
dend1 <- as.dendrogram(hcluster)
heatmap.2(netTab.seed[,SigVarsSeed.Soil], scale = "none", col = rev(bluered(100)), cexRow =0.8, cexCol = 0.8,
          Colv = dend1, Rowv=F,
          trace = "none", density.info = "none",
)

SigVarsInf.Soil <- SigVarsInf[!grepl("ASV|Stri", SigVarsInf)]

distance= dist(t(netTab.seed[,SigVarsInf.Soil]), method ="euclidian")    
hcluster = hclust(distance, method ="complete")
dend1 <- as.dendrogram(hcluster)
x<-as.matrix(netTab[,SigVarsInf.Soil, drop=FALSE])
heatmap.2(x, scale = "none", col = rev(bluered(100)), cexRow =0.8, cexCol = 0.8, Rowv=F,
          trace = "none", density.info = "none")
heatmap(x, scale = "none", col = rev(bluered(100)), Rowv = NA, Colv = NA)





CorrTab <- gjam:::.cov2Cor(outputSeed.clrNew$parameters$ematrix)
CorrTab2 <- gjam:::.cov2Cor(outputSeed.clrNew$parameters$corMu)

# parameters$corMu:This typically represents the correlation matrix between the latent factors in the model. 
# In GJAM, these latent factors are used to capture underlying patterns or structures in the data. 
# The correlation matrix provides information about how these latent factors are correlated with each other.
# parameters$ematrix: The error matrix (ematrix) usually contains information about the variances and covariances of the residuals 
# (errors) associated with the observed variables. It doesn't directly represent the correlations between the independent variables.

View(CorrTab)
dim(CorrTab)
CorrTab[]
CorrTabSigVarSoil <- data.frame(CorrTab[SigVarsSeed.Soil,!(colnames(CorrTab) %in%SigVarsSeed.Soil)])
View(CorrTabSigVarSoil)
CorrTabSigVarSoil2 <- data.frame(SoilVars = rownames(CorrTabSigVarSoil),CorrTabSigVarSoil) 



data<-data.frame(TreatGroup = rownames(netTab.seed[,SigVarsSeed.Soil]),netTab.seed[,SigVarsSeed.Soil]) %>%
  gather(key = "SoilVars", value = "RegCoeff", - TreatGroup) %>% 
  left_join(CorrTabSigVarSoil2, by = "SoilVars") %>% #View
  gather(key = "Microbes", value = "CorrCoef", -TreatGroup,-SoilVars,-RegCoeff) %>% 
  filter(abs(CorrCoef)>0.2) %>% #View
  mutate(CorrGroup = ifelse(CorrCoef >0,"Pos","Neg")) %>% 
  filter(TreatGroup == "0Seed") %>% 
  left_join(TaxInfo.df, by = c("Microbes" = "VarNames")) %>% #View
  filter(!is.na(Kingdom)) %>% 
  filter(Microbes %in% SigVarsSeed) %>% 
  group_by(TreatGroup,SoilVars,Kingdom,Phylum,Family,CorrGroup) %>% 
  summarise(Counts = n_distinct(Microbes)) %>% #View
  mutate(CorrGroup = factor(CorrGroup, levels = c("Pos", "Neg"))) %>% 
  filter(Counts>5)
  # ggplot(aes(x = interaction(CorrGroup,SoilVars),y = Counts, fill = Family,color = "black"))+
  # geom_hline(yintercept = 0, linetype = "dashed")+
  # geom_bar(position = "stack", stat = "identity")+
  # xlab("")+
  # scale_color_manual(values = "black")+
  # theme_bw()+
  # theme(panel.margin = grid::unit(-1.25, "lines"))

data <- data %>%
  mutate(interaction_var = interaction(CorrGroup, SoilVars))

ggplot(data, aes(x = interaction_var, y = Counts, fill = Family, color = CorrGroup)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("") +
  scale_color_manual(values = c("Pos" = "black", "Neg" = "white")) +  # Adjust colors as needed
  theme_bw() +
  theme(panel.margin = grid::unit(-1.25, "lines")) +
  ggtitle("Family level error covariance with Significant SEEDBANK Soil Parameters")

