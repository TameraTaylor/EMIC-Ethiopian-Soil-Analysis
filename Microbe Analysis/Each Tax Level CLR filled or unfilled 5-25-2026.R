#GLOBAL
library(tidyr)
library(ggplot2)
library(devtools)
library(dplyr)
library(compositions)
library("Tjazi")
setwd("C:/Users/small/Documents/Brady lab/Soil Paper/Luisa ASVs")


##############

#FUNGUS

##############
#1. Read in the fungus data of original taxonomy with unknowns to be just grouped together, or the created taxonomy dataset where unknowns were filled with higher known taxonomy
#1a. Read in and separate taxonomy for GROUPED UNKNOWNS
FunSeq<-read.csv("ITS-luisa-full-ave-rep.csv", header=TRUE, row.names="ASV_ID")
FUNGUS <- FunSeq %>%
  separate(taxonomy, into= c("kingdom", "phylum", "class", "order", "family", "genus"), ";")
FUNGUS<-FUNGUS[-c(49)] #remove seq
unique(FUNGUS$family) #243 unique families, 314 with unknowns
unique(FUNGUS$genus) #513 unique genus,  676 with unknowns

#FUNGUS[][FUNGUS[] == "NA"] <- NA #read in a strings, not real NAs
FUNGUS$kingdom<-gsub("k__","",as.character(FUNGUS$kingdom))
FUNGUS$phylum<-gsub("p__","",as.character(FUNGUS$phylum))
FUNGUS$class<-gsub("c__","",as.character(FUNGUS$class))
FUNGUS$order<-gsub("o__","",as.character(FUNGUS$order))
FUNGUS$family<-gsub("f__","",as.character(FUNGUS$family))
FUNGUS$genus<-gsub("g__","",as.character(FUNGUS$genus))



#1b. OR! Read in the filled unknowns dataset made in Correlation Analaysis UPDATED.R
FUNGUS<-read.csv("FUNGUS_ASVS_TaxFilled.csv", header=TRUE, row.names = "ASV_ID")



#2. AGGREGATING
Fun_phy <- FUNGUS %>%
  group_by(phylum) %>%
  summarise(across(starts_with("E"), sum, na.rm = TRUE)) %>% #sum across Soil Samples
  as.data.frame()
Fun_class <- FUNGUS %>%
  group_by(class) %>%
  summarise(across(starts_with("E"), sum, na.rm = TRUE)) %>% #sum across Soil Samples
  as.data.frame()
Fun_order <- FUNGUS %>%
  group_by(order) %>%
  summarise(across(starts_with("E"), sum, na.rm = TRUE)) %>% #sum across Soil Samples
  as.data.frame()
Fun_fam <- FUNGUS %>%
  group_by(family) %>%
  summarise(across(starts_with("E"), sum, na.rm = TRUE)) %>% #sum across Soil Samples
  as.data.frame()
Fun_gen <- FUNGUS %>%
  group_by(genus) %>%
  summarise(across(starts_with("E"), sum, na.rm = TRUE)) %>% #sum across Soil Samples
  as.data.frame()



#3. CLR TRANSFORMATION
setwd("C:/Users/small/Documents/Brady lab/Soil Paper/Luisa ASVs/Taxa level Correlations UPDATED")

CLRFun<-clr_lite(FUNGUS[1:48], samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRFun, file="Imputed CLR Unknowns Grouped/CLR_Fungus.csv")


row.names(Fun_phy)<-Fun_phy$phylum
Fun_phy<-Fun_phy[,-1]
CLRFun_phy<-clr_lite(Fun_phy, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRFun_phy, file="Imputed CLR filled Unknowns/CLR_Fungus_phylum.csv")

CLRFun_phy2<-clr_lite(Fun_phy, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRFun_phy2, file="Imputed CLR Unknowns Grouped/CLR_Fungus_phylum.csv")


row.names(Fun_class)<-Fun_class$class
Fun_class<-Fun_class[,-1]
CLRFun_class<-clr_lite(Fun_class, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRFun_class, file="Imputed CLR filled Unknowns/CLR_Fungus_class.csv")

CLRFun_class2<-clr_lite(Fun_class, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRFun_class2, file="Imputed CLR Unknowns Grouped/CLR_Fungus_class.csv")


row.names(Fun_order)<-Fun_order$order
Fun_order<-Fun_order[,-1]
CLRFun_order<-clr_lite(Fun_order, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRFun_order, file="Imputed CLR filled Unknowns/CLR_Fungus_order.csv")

CLRFun_order2<-clr_lite(Fun_order, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRFun_order2, file="Imputed CLR Unknowns Grouped/CLR_Fungus_order.csv")


row.names(Fun_fam)<-Fun_fam$family
Fun_fam<-Fun_fam[,-1]
CLRFun_fam<-clr_lite(Fun_fam, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRFun_fam, file="Imputed CLR filled Unknowns/CLR_Fungus_family.csv")

CLRFun_fam2<-clr_lite(Fun_fam, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRFun_fam2, file="Imputed CLR Unknowns Grouped/CLR_Fungus_family.csv")


row.names(Fun_gen)<-Fun_gen$genus
Fun_gen<-Fun_gen[,-1]
CLRFun_gen<-clr_lite(Fun_gen, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRFun_gen, file="Imputed CLR filled Unknowns/CLR_Fungus_genera.csv")

CLRFun_gen2<-clr_lite(Fun_gen, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRFun_gen2, file="Imputed CLR Unknowns Grouped/CLR_Fungus_genera.csv")



#####################

#BACTERIA

######################

#1a. Read in and separate taxonomy for GROUPED UNKNOWNS
BactSeq<-read.csv("C:/Users/small/Documents/Brady lab/Soil Paper/Luisa ASVs/16S-luisa-full-ave-rep.csv", header=TRUE, row.names = "ASV_ID")
BACTERIA <- BactSeq %>%
  separate(taxonomy, into= c("kingdom", "phylum", "class", "order", "family", "genus"), ";")
BACTERIA<-BACTERIA[-c(49)] #remove seq
str(BACTERIA)

#1b. OR! Read in the filled unknowns dataset made in Correlation Analaysis UPDATED.R
BACTERIA<-read.csv("C:/Users/small/Documents/Brady lab/Soil Paper/Luisa ASVs/Taxa level Correlations UPDATED/BACTERIA_ASVS_TaxFilled.csv", header=TRUE, row.names = "ASV_ID")


#2. AGGREGATING
Bact_phy <- BACTERIA %>%
  group_by(phylum) %>%
  summarise(across(starts_with("E"), sum, na.rm = TRUE)) %>% #sum across Soil Samples
  as.data.frame()
Bact_class <- BACTERIA %>%
  group_by(class) %>%
  summarise(across(starts_with("E"), sum, na.rm = TRUE)) %>% #sum across Soil Samples
  as.data.frame()
Bact_order <- BACTERIA %>%
  group_by(order) %>%
  summarise(across(starts_with("E"), sum, na.rm = TRUE)) %>% #sum across Soil Samples
  as.data.frame()
Bact_fam <- BACTERIA %>%
  group_by(family) %>%
  summarise(across(starts_with("E"), sum, na.rm = TRUE)) %>% #sum across Soil Samples
  as.data.frame()
Bact_gen <- BACTERIA %>%
  group_by(genus) %>%
  summarise(across(starts_with("E"), sum, na.rm = TRUE)) %>% #sum across Soil Samples
  as.data.frame()

#3. CLR TRANSFORMATION

CLRBact<-clr_lite(BACTERIA[1:48], samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRBact, file="Imputed CLR Unknowns Grouped/CLR_Bacteria.csv")


row.names(Bact_phy)<-Bact_phy$phylum
Bact_phy<-Bact_phy[,-1]
CLRBact_phy<-clr_lite(Bact_phy, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRBact_phy, file="Taxa level Correlations UPDATED/Imputed CLR filled Unknowns/CLR_Bacteria_phylum_.csv")

CLRBact_phy2<-clr_lite(Bact_phy, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRBact_phy2, file="Imputed CLR Unknowns Grouped/CLR_Bacteria_phylum.csv")


row.names(Bact_class)<-Bact_class$class
Bact_class<-Bact_class[,-1]
CLRBact_class<-clr_lite(Bact_class, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRBact_class, file="Taxa level Correlations UPDATED/Imputed CLR filled Unknowns/CLR_Bacteria_class_.csv")

CLRBact_class2<-clr_lite(Bact_class, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRBact_class2, file="Imputed CLR Unknowns Grouped/CLR_Bacteria_class.csv")


row.names(Bact_order)<-Bact_order$order
Bact_order<-Bact_order[,-1]
CLRBact_order<-clr_lite(Bact_order, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRBact_order, file="Taxa level Correlations UPDATED/Imputed CLR filled Unknowns/CLR_Bacteria_order_.csv")

CLRBact_order2<-clr_lite(Bact_order, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRBact_order2, file="Imputed CLR Unknowns Grouped/CLR_Bacteria_order.csv")


row.names(Bact_fam)<-Bact_fam$family
Bact_fam<-Bact_fam[,-1]
CLRBact_fam<-clr_lite(Bact_fam, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRBact_fam, file="Taxa level Correlations UPDATED/Imputed CLR filled Unknowns/CLR_Bacteria_family_.csv")

CLRBact_fam2<-clr_lite(Bact_fam, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRBact_fam2, file="Imputed CLR Unknowns Grouped/CLR_Bacteria_family.csv")


row.names(Bact_gen)<-Bact_gen$genus
Bact_gen<-Bact_gen[,-1]
CLRBact_gen<-clr_lite(Bact_gen, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRBact_gen, file="Taxa level Correlations UPDATED/Imputed CLR filled Unknowns/CLR_Bacteria_genera.csv")

CLRBact_gen2<-clr_lite(Bact_gen, samples_are = "cols", method="logunif", replicates=1000)
write.csv(CLRBact_gen2, file="Imputed CLR Unknowns Grouped/CLR_Bacteria_genera.csv")
