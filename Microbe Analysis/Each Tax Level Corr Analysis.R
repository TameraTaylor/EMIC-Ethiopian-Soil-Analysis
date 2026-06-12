#GLOBAL
library(tidyr)
library(ggplot2)
library(devtools)
library(dplyr)
library("psych")
library(corrplot)
library(pheatmap)
library("openxlsx")
setwd("C:/Users/small/Documents/Brady lab/Soil Paper/Taxa Level Correlations UPDATED")

######################
##Striga Data
#######################

#Get Striga data and add to datasets
d<-read.csv("C:/Users/small/Documents/Brady lab/Soil Paper/FullSoilV7.csv", row.names = "ID")
d<-d[,c(38:45)] #only striga columns
d<-d[-c(49, 50),]#remove e49 and 50 (no microbe data) 
data_norm<-d[,c(5:8)]
data_raw<-d[,c(1:4)]




####################
#Create Functions for adding Striga data and running correlation
####################

run_corr_analysis <- function(microbe_data, out_path) {
  data_matrix <- rbind(microbe_data, t(data_norm))
  dataT <- t(data_matrix) %>% as.matrix()
  n_traits <- 4
  tax_cols <- 1:(ncol(dataT) - n_traits)
  trait_cols <- (ncol(dataT) - n_traits + 1):ncol(dataT)

  cor2 <- corr.test(dataT[, tax_cols], dataT[, trait_cols],
                    method = "pearson", adjust = "bonferroni")

  cor2sigr <- ifelse(cor2$p < 0.05, cor2$r, NA)
  cor2sigrNArem <- cor2sigr[apply(!is.na(cor2sigr), 1, any), , drop = FALSE]

  wb <- createWorkbook()
  addWorksheet(wb, "P")
  addWorksheet(wb, "Padj")
  addWorksheet(wb, "R")
  addWorksheet(wb, "SigR-P")
  writeData(wb, "P", cor2$p, rowNames=TRUE)
  writeData(wb, "Padj", cor2$p.adj, rowNames=TRUE)
  writeData(wb, "R", cor2$r, rowNames=TRUE)
  writeData(wb, "SigR-P", cor2sigrNArem, rowNames=TRUE)
  saveWorkbook(wb, out_path, overwrite = TRUE)
}


####################
#Read in Microbe datas
####################






########################
#Run Correlations
########################
setwd("C:/Users/small/Documents/Brady lab/Soil Paper/Luisa ASVs/Taxa level Correlations UPDATED")

run_corr_analysis(CLRFun_gen2, "Imputed CLR Unknowns grouped/CORR_Fungus_genera.xlsx")
run_corr_analysis(CLRFun_fam2, "Imputed CLR Unknowns grouped/CORR_Fungus_family.xlsx")
run_corr_analysis(CLRFun_order2, "Imputed CLR Unknowns grouped/CORR_Fungus_order.xlsx")
run_corr_analysis(CLRFun_class2, "Imputed CLR Unknowns grouped/CORR_Fungus_class.xlsx")
run_corr_analysis(CLRFun_phy2, "Imputed CLR Unknowns grouped/CORR_Fungus_phylum.xlsx")

run_corr_analysis(CLRFun_gen, "Unknowns filled/CORR_Fungus_genera_filledUnknowns.xlsx")
run_corr_analysis(CLRFun_fam, "Unknowns filled/CORR_Fungus_family_filledUnknowns.xlsx")
run_corr_analysis(CLRFun_order, "Unknowns filled/CORR_Fungus_order_filledUnknowns.xlsx")
run_corr_analysis(CLRFun_class, "Unknowns filled/CORR_Fungus_class_filledUnknowns.xlsx")
run_corr_analysis(CLRFun_phy, "Unknowns filled/CORR_Fungus_genera_filledUnknowns.xlsx")

run_corr_analysis(CLRBact_gen2, "Imputed CLR Unknowns grouped/CORR_Bacteria_genera.xlsx")
run_corr_analysis(CLRBact_fam2, "Imputed CLR Unknowns grouped/CORR_Bacteria_family.xlsx")
run_corr_analysis(CLRBact_order2, "Imputed CLR Unknowns grouped/CORR_Bacteria_order.xlsx")
run_corr_analysis(CLRBact_class2, "Imputed CLR Unknowns grouped/CORR_Bacteria_class.xlsx")
run_corr_analysis(CLRBact_phy2, "Imputed CLR Unknowns grouped/CORR_Bacteria_phylum.xlsx")

run_corr_analysis(CLRBact_gen, "Unknowns filled/CORR_Bacteria_genera_filledUnknowns.xlsx")
run_corr_analysis(CLRBact_fam, "Unknowns filled/CORR_Bacteria_family_filledUnknowns.xlsx")
run_corr_analysis(CLRBact_order, "Unknowns filled/CORR_Bacteria_order_filledUnknowns.xlsx")
run_corr_analysis(CLRBact_class, "Unknowns filled/CORR_Bacteria_class_filledUnknowns.xlsx")
run_corr_analysis(CLRBact_phy, "Unknowns filled/CORR_Bacteria_genera_filledUnknowns.xlsx")






##################
#Order the genera for better figure making
#################

setwd("C:/Users/small/Documents/Brady lab/Soil Paper/Luisa ASVs/Taxa level Correlaions UPDATED")


#read in cor2sigr and its matching original tax file that has larger tax info
cor2sigrNArem<-read.xlsx("Imputed CLR Unknowns grouped/CORR_Fungus_genera.xlsx", "SigR-P", rowNames = TRUE)
data<-read.csv("FUNGUS ASVs.csv", row.names="ASV_ID")

cor2sigrNArem<-read.xlsx("Unknowns filled/CORR_Fungus_genera_filledUnknowns.xlsx", "SigR-P", rowNames = TRUE)
data<-read.csv("FUNGUS_ASVS_TaxFilled.csv", row.names = "X")

cor2sigrNArem<-read.xlsx("Imputed CLR Unknowns grouped/CORR_Bacteria_genera.xlsx", "SigR-P", rowNames = TRUE)
data<-read.csv("BACTERIA ASVs.csv", row.names = "ASV_ID")

cor2sigrNArem<-read.xlsx("Unknowns filled/CORR_Bacteria_genera_filledUnknowns.xlsx", "SigR-P", rowNames = TRUE)
data<-read.csv("BACTERIA_ASVS_TaxFilled.csv", row.names = "X")


df<-row.names(cor2sigrNArem) 


#Find matches and reorder
data_filtered <- data[data$genus %in% df, ]
data_filtered<-data_filtered[with(data_filtered, order(phylum, class, order, family)), ] #order with multiple nested layers
ordered<-data_filtered$genus %>% unique() #get list of unique genera in their nested order

cor2sigrNArem2 <- cor2sigrNArem[match(ordered, row.names(cor2sigrNArem)), ] #reorder the genera based on ordered list
#cor2sigrNArem2 <-rbind(cor2sigrNArem2, cor2sigrNArem[31,])
str(cor2sigrNArem)#check str of cor2sigrNArem og and 2, should have same # of rows
str(cor2sigrNArem2)

#graph
pheatmap(cor2sigrNArem2, fontsize_row = 6, fontsize_col = 6, cluster_cols = FALSE, cluster_rows = FALSE, 
         na_col = "grey", color=colorRampPalette(c("red", "white", "blue"))(50), main = "Bacterial Genera")

data_filtered<-data_filtered %>% distinct(genus, .keep_all = TRUE)
write.csv(data_filtered[c(49:54)], "Imputed CLR Unknowns grouped/SigFungusTaxOrdered.csv")
