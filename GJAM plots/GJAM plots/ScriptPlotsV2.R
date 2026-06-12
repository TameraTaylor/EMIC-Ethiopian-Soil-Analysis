#Files to share with Tamera so she can have the number she want

save(pS.seed, fungS.seed, SoilS.seed, StrigaInf.seed,netTab.seed,SigVars.seed,
     TreatTab.seed,Avg_Abund.seed,SignInfo.seed,Seq16S,
     divnet.shannon, divnet.shannon.Fun, # Diversity index calculated by bacteria and fungi
     xdata.Seed,Striga_seedbank, TaxInfoBac,SigVars,
     TreatTab,Avg_Abund,SignInfo,TaxInfo.df,
     SensTabInf,pS, fungS, SoilS, netTab, SigVarsInf, #GJAM model results using striga infestation as the explanatory variable
     xdata,Suic_germ, #Data from Suicidal germination as the explanatory variable
     file = "GJAMResults_ML.RData")

#Script used for generating the plots

# impact of diversity


#Bacteria plots

data.frame(xdata,divnet.bact.shannon,Suic_germ) %>% 
  filter(Sample.code != "E15") %>% 
  filter(Sample.code != "E48") %>% 
  ggplot(aes(x = Suicidal.germination.....Mean, y = estimate, col = GroupVar))+
  geom_errorbarh(aes(xmax = Suicidal.germination.....Mean+Suicidal.germination.....SE,
                     xmin = Suicidal.germination.....Mean - Suicidal.germination.....SE, height = .02))+
  ylab("Estimate of Shannon index (Bacteria)")+
  xlab("Suicidal germination (%) Mean")+
  geom_point()+
  theme_bw()

data.frame(xdata.Inf,divnet.bact.shannon,Striga_field_inf) %>% 
  filter(Sample.code != "E15") %>% 
  filter(Sample.code != "E48") %>% 
  ggplot(aes(x = Striga.count.per.m2.Mean, y = estimate, col = GroupVar))+
  geom_errorbarh(aes(xmax = Striga.count.per.m2.Mean+Striga.count.per.m2.SE,
                     xmin = Striga.count.per.m2.Mean - Striga.count.per.m2.SE, height = .02))+
  ylab("Estimate of Shannon index (Bacteria)")+
  xlab("Striga count per m²")+
  geom_point()+
  theme_bw()

data.frame(xdata.Inf,divnet.bact.shannon,Striga_seedbank) %>% 
  filter(Sample.code != "E15") %>% 
  filter(Sample.code != "E48") %>% 
  ggplot(aes(x = Striga.count.per.m2.Mean, y = estimate, col = GroupVar))+
  geom_errorbarh(aes(xmax = Striga.count.per.m2.Mean+Striga.count.per.m2.SE,
                     xmin = Striga.count.per.m2.Mean - Striga.count.per.m2.SE, height = .02))+
  ylab("Estimate of Shannon index (Bacteria)")+
  xlab("Striga count per m²")+
  geom_point()+
  theme_bw()


#Fungi

data.frame(xdata,divnet.shannon.Fun,Suic_germ) %>% 
  filter(Sample.code != "E15") %>% 
  filter(Sample.code != "E48") %>% 
  ggplot(aes(x = Suicidal.germination.....Mean, y = estimate, col = GroupVar))+
  geom_errorbarh(aes(xmax = Suicidal.germination.....Mean+Suicidal.germination.....SE,
                     xmin = Suicidal.germination.....Mean - Suicidal.germination.....SE, height = .02))+
  ylab("Estimate of Shannon index (Fungi)")+
  xlab("Suicidal germination (%) Mean")+
  geom_point()+
  theme_bw()

data.frame(xdata.Inf,divnet.shannon.Fun,Striga_field_inf) %>% 
  filter(Sample.code != "E15") %>% 
  filter(Sample.code != "E48") %>% 
  ggplot(aes(x = Striga.count.per.m2.Mean, y = estimate, col = GroupVar))+
  geom_errorbarh(aes(xmax = Striga.count.per.m2.Mean+Striga.count.per.m2.SE,
                     xmin = Striga.count.per.m2.Mean - Striga.count.per.m2.SE, height = .02))+
  ylab("Estimate of Shannon index (Fungi)")+
  xlab("Striga count per m²")+
  # geom_smooth(aes(group = 1))+
  geom_point()+
  theme_bw()


# Analysis GJAM infestation

#Sensitivity
SensTabInf <- data.frame(GroupVar = c(rep("proc", nrow(pS)), rep("fung", nrow(fungS)), rep("Soil factors", nrow(SoilS)),
                                      rep("Striga", nrow(StrigaSeedS))),
                         rbind(pS, fungS, SoilS, StrigaSeedS))
View(SensTabInf)
SensTabInf %>% 
  gather(key = "TreatVar", value = "Sensitivity", -GroupVar) %>% 
  mutate(GroupVar = factor(GroupVar, levels = c("proc", "fung", "Soil factors","Striga"))) %>% 
  mutate(TreatVar = gsub("GroupVar", "", TreatVar)) %>% 
  mutate(TreatVar = factor(TreatVar, levels = c("0Inf", "Inf10", "Inf1040", "Inf4080", "Inf80180"))) %>% #summary
  ggplot(aes(x = TreatVar, y = Sensitivity, fill = GroupVar))+
  geom_boxplot()+
  scale_y_log10() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),#legend.position = "none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
colnames(SensTab)

SensTabInf %>% 
  gather(key = "TreatVar", value = "Sensitivity", -GroupVar) %>% 
  filter(GroupVar != "Striga") %>% 
  mutate(GroupVar = factor(GroupVar, levels = c("proc", "fung", "Soil factors"))) %>% 
  mutate(TreatVar = gsub("GroupVar", "", TreatVar)) %>% 
  mutate(TreatVar = factor(TreatVar, levels = c("0Inf", "Inf10", "Inf1040", "Inf4080", "Inf80180"))) %>% #summary
  ggplot(aes(x = TreatVar, y = Sensitivity, fill = GroupVar))+
  geom_boxplot()+
  scale_y_log10() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),#legend.position = "none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#Heatmap

#Coefficients
heatmap(selected.model$parameters$betaStandXWmu)
library("gplots")
netTab <- selected.model$parameters$betaStandXWmu
rownames(netTab) <- gsub("GroupVar", "", rownames(netTab))
hT <- heatmap.2(netTab, scale = "none", col = rev(bluered(100)), cexRow =0.8,#lhei=lhei, lwid=lwid,lmat =lmat,
                hclustfun = function(x) hclust(x,method = "average"),
                trace = "none", density.info = "none"#, labCol = LabelsCode[colnames(netTab),]$AbbTax
)
dev.off()

hTInf <- heatmap.2(netTab[,SigVarsInf], scale = "none", col = rev(bluered(100)), cexRow =0.8,#lhei=lhei, lwid=lwid,lmat =lmat,
                   hclustfun = function(x) hclust(x,method = "average"),
                   trace = "none", density.info = "none"#, labCol = LabelsCode[colnames(netTab),]$AbbTax
)
dev.off()

SigVars %>% 
  left_join(TaxInfoBac)


#Bubble plot
View(SigVars)

TreatTab <- data.frame(VarNames = colnames(netTab),t(netTab))
View(TaxInfo.df)


library(ggrepel)

TreatTab %>% 
  left_join(Avg_Abund, by = "VarNames") %>% 
  mutate(SigVar = ifelse(VarNames %in% SigVarsInf, "Sig", "NS")) %>% 
  gather(key = "TreatVar", value = "Coeff", -CLR_Abund, -VarNames, -SigVar) %>%
  mutate(TreatVar = gsub("X0", "0", TreatVar)) %>% 
  left_join(SignInfo, by= c("VarNames", "TreatVar")) %>% 
  mutate(sig95 = ifelse(is.na(sig95), "ns", sig95)) %>% 
  left_join(TaxInfo.df, by = "VarNames") %>% 
  mutate(Phylum = factor(Phylum, levels = unique(Phylum))) %>%
  filter(!is.na(Kingdom)) %>% 
  group_by(Kingdom, TreatVar) %>%
  mutate(TaxLab = ifelse(Coeff %in% tail(sort(Coeff), 3), paste(Genus, Species), NA)) %>% #View
  mutate(TaxLab = ifelse(Coeff %in% head(sort(Coeff), 3), paste(Genus, Species), TaxLab)) %>% #View
  mutate(TaxLab = ifelse(TaxLab == "NA NA", Family, TaxLab)) %>% #View
  mutate(TaxLab = ifelse(TaxLab == "NA", paste(Kingdom, Phylum), TaxLab)) %>%
  # mutate(TaxLab = ifelse(TaxLab == "NA", Class, TaxLab)) %>%
  # mutate(TaxLab = ifelse(TaxLab == "NA", Phylum, TaxLab)) %>% View
  ggplot(aes(x = CLR_Abund, y = Coeff, color = Phylum, alpha = sig95))+
  geom_point()+
  geom_text_repel(aes(label = TaxLab))+
  facet_grid(Kingdom~TreatVar, scales = "free_y")+
  scale_alpha_manual(values = c(0.7,0.05))+
  theme_bw()


# Analysis GJAM seed bank

#Sensitivity

SensTabInf.seed <- data.frame(GroupVar = c(rep("proc", nrow(pS.seed)), rep("fung", nrow(fungS.seed)), 
                                           rep("Soil factors", nrow(SoilS.seed)),
                                           rep("Striga", nrow(StrigaInf.seed))),
                              rbind(pS.seed, fungS.seed, SoilS.seed, StrigaInf.seed))
View(SensTabInf.seed)
colnames(SensTabInf.seed)
SensTabInf.seed %>% 
  gather(key = "TreatVar", value = "Sensitivity", -GroupVar) %>% 
  mutate(GroupVar = factor(GroupVar, levels = c("proc", "fung", "Soil factors","Striga"))) %>% 
  mutate(TreatVar = gsub("GroupVar", "", TreatVar)) %>% 
  mutate(TreatVar = factor(TreatVar, levels = c("0Seed", "Seed10", "Seed1020", "Seed2061"))) %>% #summary
  ggplot(aes(x = TreatVar, y = Sensitivity, fill = GroupVar))+
  geom_boxplot()+
  scale_y_log10() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),#legend.position = "none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
colnames(SensTab)

SensTabInf.seed %>% 
  gather(key = "TreatVar", value = "Sensitivity", -GroupVar) %>% 
  filter(GroupVar != "Striga") %>% 
  mutate(GroupVar = factor(GroupVar, levels = c("proc", "fung", "Soil factors"))) %>% 
  mutate(TreatVar = gsub("GroupVar", "", TreatVar)) %>% 
  mutate(TreatVar = factor(TreatVar, levels = c("0Seed", "Seed10", "Seed1020", "Seed2061"))) %>% #summary
  ggplot(aes(x = TreatVar, y = Sensitivity, fill = GroupVar))+
  geom_boxplot()+
  scale_y_log10() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),#legend.position = "none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#Heatmap

#Coefficients
library("gplots")
hT.seed <- heatmap.2(netTab.seed, scale = "none", col = rev(bluered(100)), cexRow =0.8,#lhei=lhei, lwid=lwid,lmat =lmat,
                     hclustfun = function(x) hclust(x,method = "average"),
                     trace = "none", density.info = "none"#, labCol = LabelsCode[colnames(netTab),]$AbbTax
)
dev.off()

hTSigSeed <- heatmap.2(netTab.seed[,SigVarsSeed], scale = "none", col = rev(bluered(100)), cexRow =0.8,#lhei=lhei, lwid=lwid,lmat =lmat,
                       # hclustfun = function(x) hclust(x,method = "average"),
                       trace = "none", density.info = "none"#, labCol = LabelsCode[colnames(netTab),]$AbbTax
)
dev.off()


#Bubble plot

library(ggrepel)

TreatTab.seed %>% 
  left_join(Avg_Abund.seed, by = "VarNames") %>% 
  mutate(SigVar = ifelse(VarNames %in% SigVarsInf, "Sig", "NS")) %>% 
  gather(key = "TreatVar", value = "Coeff", -CLR_Abund, -VarNames, -SigVar) %>%
  mutate(TreatVar = gsub("X0", "0", TreatVar)) %>% 
  left_join(SignInfo.seed, by= c("VarNames", "TreatVar")) %>% 
  mutate(sig95 = ifelse(is.na(sig95), "ns", sig95)) %>% 
  left_join(TaxInfo.df, by = "VarNames") %>% 
  mutate(Phylum = factor(Phylum, levels = unique(Phylum))) %>%
  filter(!is.na(Kingdom)) %>% 
  group_by(Kingdom, TreatVar) %>%
  mutate(TaxLab = ifelse(Coeff %in% tail(sort(Coeff), 3), paste(Genus, Species), NA)) %>% #View
  mutate(TaxLab = ifelse(Coeff %in% head(sort(Coeff), 3), paste(Genus, Species), TaxLab)) %>% #View
  mutate(TaxLab = ifelse(TaxLab == "NA NA", Family, TaxLab)) %>% #View
  mutate(TaxLab = ifelse(TaxLab == "NA", paste(Kingdom, Phylum), TaxLab)) %>%
  # mutate(TaxLab = ifelse(TaxLab == "NA", Class, TaxLab)) %>%
  # mutate(TaxLab = ifelse(TaxLab == "NA", Phylum, TaxLab)) %>% View
  ggplot(aes(x = CLR_Abund, y = Coeff, color = Phylum, alpha = sig95))+
  geom_point()+
  geom_text_repel(aes(label = TaxLab))+
  facet_grid(Kingdom~TreatVar, scales = "free_y")+
  scale_alpha_manual(values = c(0.7,0.05))+
  theme_bw()


# Diversity related to the seedbank

data.frame(xdata.Seed,divnet.shannon,Striga_seedbank) %>% #colnames()
  ggplot(aes(x = Striga.seeds.per.150.g.of.soil.Mean, y = estimate, col = GroupVar, size = Striga.seeds.per.150.g.of.soil.SE))+
  geom_point()

##################################################################################################

# New Files as requested by Tamera

save(divnet.bact.shannon,
     Striga_field_inf,
     xdata.Inf,
     # Sample.code, #this is the Sample Code as present in the xdata
     TaxInfoBac, TaxInfoFun,
     SigVarsSeed,
     # SigVars.Suic, # could not find it easily and I think you will not need that anyway so if you really want let me know and I can generate the file again
     file = "GJAMExtraObjects_ML.RData")

#Sensitivity Analysis

# Analysis GJAM infestation

#Sensitivity analysis

selected.model <- outputInf.clr # this is an example of GJAM model that I use to plot. So I can reuse my script. This is what I mean by a selected.model
selected.model$inputs$y
ynames <- colnames(selected.model$inputs$y)
rev(ynames)
ynames
View(TaxInfo[TaxInfo$TaxCode %in% ynames,])

proc  <- ynames[1:882]
fung  <- ynames[883:1450]
Soil  <- ynames[1451:1480]
StrigaSeed  <- ynames[1481:1482]


pS   <- gjamSensitivity(selected.model, proc)
fungS <- gjamSensitivity(selected.model, fung)
SoilS   <- gjamSensitivity(selected.model, Soil)
StrigaSeedS   <- gjamSensitivity(selected.model, StrigaSeed)
save.image("AnalysisGetahuni7.RData")


SensTabInf <- data.frame(GroupVar = c(rep("proc", nrow(pS)), rep("fung", nrow(fungS)), rep("Soil factors", nrow(SoilS)),
                                      rep("Striga", nrow(StrigaSeedS))),
                         rbind(pS, fungS, SoilS, StrigaSeedS))
View(SensTabInf)
SensTabInf %>% 
  gather(key = "TreatVar", value = "Sensitivity", -GroupVar) %>% 
  mutate(GroupVar = factor(GroupVar, levels = c("proc", "fung", "Soil factors","Striga"))) %>% 
  mutate(TreatVar = gsub("GroupVar", "", TreatVar)) %>% 
  mutate(TreatVar = factor(TreatVar, levels = c("0Inf", "Inf10", "Inf1040", "Inf4080", "Inf80180"))) %>% #summary
  ggplot(aes(x = TreatVar, y = Sensitivity, fill = GroupVar))+
  geom_boxplot()+
  scale_y_log10() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),#legend.position = "none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
colnames(SensTab)

SensTabInf %>% 
  gather(key = "TreatVar", value = "Sensitivity", -GroupVar) %>% 
  filter(GroupVar != "Striga") %>% 
  mutate(GroupVar = factor(GroupVar, levels = c("proc", "fung", "Soil factors"))) %>% 
  mutate(TreatVar = gsub("GroupVar", "", TreatVar)) %>% 
  mutate(TreatVar = factor(TreatVar, levels = c("0Inf", "Inf10", "Inf1040", "Inf4080", "Inf80180"))) %>% #summary
  ggplot(aes(x = TreatVar, y = Sensitivity, fill = GroupVar))+
  geom_boxplot()+
  scale_y_log10() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),#legend.position = "none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Making TreatTab
netTab <- selected.model$parameters$betaStandXWmu
rownames(netTab) <- gsub("GroupVar", "", rownames(netTab))
TreatTab <- data.frame(VarNames = colnames(netTab),t(netTab))


# GJAM for the Striga Seedbank - This is how I build the objects and run the model


ydataSeed.clr <- cbind(clrdata.list$Bacteria, clrdata.list$Fungi,
                       as.matrix(scale(Soil_Info[,-1])),#scale(Striga_field_inf[,6:7]),
                       scale(Striga_field_inf[,6:7]))

dim(ydataSeed.clr);dim(xdata.Seed)
mlSeed.clr <- list( ng = 10000, burnin = 2000, typeNames = rep("CON",ncol(ydataSeed.clr)), 
                    # reductList = list(r = 8, N = 15),
                    random = "District" )
outputSeed.clr <- gjam(~ GroupVar, data.frame(xdata.Seed), data.frame(ydataSeed.clr), 
                       modelList = ml.clr)

# Then whenever I am runing a script a select a model which is why I have a selected.model. So I can reuse my code.