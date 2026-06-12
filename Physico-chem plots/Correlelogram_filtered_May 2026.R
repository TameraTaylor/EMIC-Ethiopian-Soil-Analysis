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
#PEARSON
data_ready <- as.matrix(data_ready)
C <- cor(data_ready,  use="pairwise.complete.obs") #uses Pearson by default
cortest<-corr.test(data_ready, adjust= "bonferroni") #default is Pearson
cortestp<-cortest$p
cortestr<-cortest$r #its important to check that C and cortestr are identical!
write.csv(C, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsPV3cor.csv")
write.csv(cortestp, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsPV3corrtestPvalues.csv")
write.csv(cortestr, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsPV3corrtestRvalues.csv")

#SPEARMAN
C <- cor(data_ready,  use="pairwise.complete.obs", method = "spearman") 
cortest<-corr.test(data_ready, method = "spearman", adjust= "bonferroni") 
cortestp<-cortest$p
cortestr<-cortest$r #its important to check that C and cortestr are identical!
write.csv(C, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsSV3cor.csv")
write.csv(cortestp, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsSV3corrtestPvalues.csv")
write.csv(cortestr, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsSV3corrtestRvalues.csv")

#6-6-2026
#Should transform the compositional parts of the data first, then run correlations
library(compositions)
#SOIL TEXTURE
#normalize
texture<-c("Clay", "Silt", "Sand", "organic", "Carbonated_lime")
data_ready[, texture]<- data_ready[, texture] %>%
  apply(1, function(x) x / sum(x) * 100) %>%
  t()
#transform
comp_subset <- data_ready[, texture] 
clr_transformed <- as.data.frame(clr(comp_subset+0.001))
#replace
data_ready[, texture] <- clr_transformed

#NUTRIENT PROFILE
#normalize
nutrient<-c("Mg", "Na", "K", "Ca", "C")
data_ready[, nutrient]<- data_ready[, nutrient] %>%
  apply(1, function(x) x / sum(x) * 100) %>%
  t()
#transform
comp_subset <- data_ready[, nutrient] 
clr_transformed <- as.data.frame(clr(comp_subset+0.001))
#replace
data_ready[, nutrient] <- clr_transformed

###Then run Correlations
#SPEARMAN
C <- cor(data_ready,  use="pairwise.complete.obs", method = "spearman") 
cortest<-corr.test(data_ready, method = "spearman", adjust= "bonferroni") 
cortestp<-cortest$p
cortestr<-cortest$r #its important to check that C and cortestr are identical!
write.csv(C, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsSV4cor.csv")
write.csv(cortestp, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsSV4corrtestPvalues.csv")
write.csv(cortestr, "C:/Users/small/Documents/Brady lab/Soil Paper/Correlogram/CorrelationsSV4corrtestRvalues.csv")



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
         pch.cex=.8,
         tl.cex = .5)
dev.off()




