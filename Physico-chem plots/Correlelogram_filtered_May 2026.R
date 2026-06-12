#Remove suicidal germ, site, altitude, lon & lat
#remove Striga_count and Sorghum_counts, onlye used Normalized Striga Count
#Reorder axis so striga stuff is together
library(dplyr)
library(janitor)
library("psych") 
library(corrplot)

###
# read in data
###

setwd("C:/Users/small/Documents/Brady lab/Soil Paper")
data<-read.csv("FullSoilV3.csv")

data_no_extra <- data %>% select(!c(site:lat,Striga_Count_Ave:Striga_count.Sorghum_count,
                                    Suicidal_Germ_Ave:Suicidal_Germ_SE))



###
# make correlation
###

data_no_extra[is.na(data_no_extra)] = 0 #change NAs to 0s.
data_no_extra<-remove_constant(data_no_extra, na.rm = FALSE, quiet = TRUE)
rownames(data_no_extra)<-data_no_extra[,1]
data_no_extra<-data_no_extra[,-1]

 

data_ready <- data_no_extra %>% select("pH","Mg_avail","Mg_total","Mg","Na_avail","Na_total","Na","K_avail","K_total","K","Ca_avail"
                                       ,"Ca_total","Ca","P_avail","P_total","S_avail","S_total","S.delivery","C.S_ratio","C.OS.ratio","C"
                                       ,"C_inorganic","N","C.N_ratio","N.delivery","Microbial_activity","organic","Carbonated_lime","Clay"
                                       ,"Silt","Sand","Clay_humus","Observed","Chao1","se.chao1","ACE","se.ACE","Shannon","Simpson","InvSimpson"
                                       ,"Fisher","Abditibacteriota","Acidobacteriota","Actinobacteriota","Armatimonadota","Bacteroidota","Bdellovibrionota"
                                       ,"Chloroflexi","Cyanobacteria","Deinococcota","Entotheonellaeota","Fibrobacterota","Firmicutes","Gemmatimonadota"
                                       ,"Myxococcota","Nitrospirota","Patescibacteria","Planctomycetota","Proteobacteria","Spirochaetota"
                                       ,"Sumerlaeota","Verrucomicrobiota","Normalized_Striga_Count_Ave","Normalized_Striga_Count_SE"
                                       ,"Seedbank_Ave","Seedbank_SE")
#other version without microbe
data_ready <- data_no_extra %>% select("pH","Mg_avail","Mg_total","Mg","Na_avail","Na_total","Na","K_avail","K_total","K","Ca_avail"
                                       ,"Ca_total","Ca","P_avail","P_total","S_avail","S_total","S.delivery","C.S_ratio","C.OS.ratio","C"
                                       ,"C_inorganic","N","C.N_ratio","N.delivery","Microbial_activity","organic","Carbonated_lime","Clay"
                                       ,"Silt","Sand","Clay_humus","Normalized_Striga_Count_Ave","Normalized_Striga_Count_SE"
                                       ,"Seedbank_Ave","Seedbank_SE")

#V3 third version with shannon diversity
setwd("C:/Users/small/Documents/Brady lab/Soil Paper")
data<-read.csv("FullSoilV3.3.csv")
data_no_extra <- data %>% select(!c(site:lat,Striga_Count_Ave:Striga_count.Sorghum_count))
data_ready <- data_no_extra %>% select("pH","Mg_avail","Mg_total","Mg","Na_avail","Na_total","Na","K_avail","K_total","K","Ca_avail"
                                       ,"Ca_total","Ca","P_avail","P_total","S_avail","S_total","S.delivery","C.S_ratio","C.OS.ratio","C"
                                       ,"C_inorganic","N","C.N_ratio","N.delivery","Microbial_activity","organic","Carbonated_lime","Clay"
                                       ,"Silt","Sand","Clay_humus","BACTERIA_ShannonIndex", "FUNGUS_ShannonIndex"
                                       ,"log.Infestation_Ave",	"log.Infestation_SE", "log.Seedbank_Ave","log.Seedbank_SE",	)

#9-27-24 Need to remove 15 and 48! If skip the NA to 0 part (lines 26 & 27), should be removed from cor. 
#Looks like it doesn't affect the actual figure, just the numbers. 

data_ready <- as.matrix(data_ready)
C <- cor(data_ready,  use="pairwise.complete.obs") #uses Pearson by default
cortest<-corr.test(data_ready, adjust= "bonferroni") #default is Pearson
cortestp<-cortest$p
cortestr<-cortest$r #its important to check that C and cortestr are identical!
write.csv(C, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsPV3cor.csv")
write.csv(cortestp, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsPV3corrtestPvalues.csv")
write.csv(cortestr, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsPV3corrtestRvalues.csv")

C <- cor(data_ready,  use="pairwise.complete.obs", method = "spearman") 
cortest<-corr.test(data_ready, method = "spearman", adjust= "bonferroni") 
cortestp<-cortest$p
cortestr<-cortest$r #its important to check that C and cortestr are identical!
write.csv(C, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsSV3cor.csv")
write.csv(cortestp, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsSV3corrtestPvalues.csv")
write.csv(cortestr, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsSV3corrtestRvalues.csv")


###
# make pretty plots
###
dev.new()
corrplot(C,type='lower',  
         method = 'color', 
         diag = FALSE, 
         p.mat = cortestp, 
         insig = 'label_sig', 
         sig.level = c(0.001, 0.01, 0.05), 
         order="original", 
         tl.col="black",
         pch.cex=1,
         tl.cex = 1)
dev.off()



#######################
#5-5-2025
#Using most updated soil data v10 with many taxa
######################
setwd("C:/Users/small/Documents/Brady lab/Soil Paper")
data<-read.csv("FullSoilV10.csv")
data_no_extra <- data %>% select(!c(site:lat,Normalized_Striga_Count_Ave:Seedbank_SE,Agro_Climate_type,
                                    Abditibacteriota:Unknown_Fungus))
rownames(data_no_extra)<-data_no_extra[,1]
data_no_extra<-data_no_extra[,-1]


data_ready <- as.matrix(data_no_extra)
C <- cor(data_ready,  use="pairwise.complete.obs") #uses Pearson by default
cortest<-corr.test(data_ready, method="pearson", adjust="bonferroni") 
cortestp<-cortest$p
cortestr<-cortest$r #its important to check that C and cortestr are identical!
cortestpadj<-cortest$p.adj

library("openxlsx")
# Create a new workbook and add a sheet
wb <- createWorkbook()
addWorksheet(wb, "P")
addWorksheet(wb, "R")
# Write data to the sheet
writeData(wb, "P", cortestp, rowNames=TRUE)
writeData(wb, "R", cortestr, rowNames=TRUE)
# Save the workbook
saveWorkbook(wb, "Correlogram/CorStats.xlsx", overwrite = TRUE)

dev.new()
g <- corrplot(C,type='lower',  
              method = 'color', 
              diag = FALSE, 
              p.mat = cortestp, 
              insig = 'label_sig', 
              sig.level = c(0.001, 0.01, 0.05), 
              order="original", 
              tl.col="black",
              pch.cex=1,
              tl.cex = .5)
dev.off()
