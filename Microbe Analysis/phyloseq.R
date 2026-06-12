######################
# Transfer data to Phyloseq #
#########################

setwd("C:/Users/small/Documents/Brady lab/Soil Paper/Luisa ASVs")
library(tidyr)
library(ggplot2)
library(dplyr)
library("phyloseq")
library(vegan)
#-----

soil<-read.csv("FullSoilV5.csv", row.names = "ID")
soil<-soil[c(1:4)]


#FUNGUS######################
data<-read.csv("Taxa level Correlations UPDATED/FUNGUS ASVs.csv", row.names = "ASV_ID") #Fungus reps averaged
otumat<-data[1:48]
otumat[] <-lapply(otumat, as.integer)
otumat<-as.matrix(otumat)

taxmat<-data[49:54]
taxmat<-as.matrix(taxmat)

#make into phyloseq object
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SAM=sample_data(soil) #this one can be a dataframe
physeq = phyloseq(OTU, TAX, SAM)
physeq
 
#Rarefaction
sample_sums(physeq)
#rarecurve(t(otumat), step=50, cex=0.5)
ps.rarefied = rarefy_even_depth(physeq, rngseed=1, sample.size=0.9*min(sample_sums(physeq)), replace=F)
sample_sums(ps.rarefied)
write.csv(otu_table(ps.rarefied), "Fungus rarefied counts.csv")
plot_bar(physeq, fill = "phylum")
plot_bar(ps.rarefied, fill = "phylum")
otu_matrix <- as(otu_table(physeq), "matrix")
rarecurve(t(otu_matrix), step = 100, label = T)

diversity1<-estimate_richness(physeq)
diversity2<-estimate_richness(ps.rarefied)
write.csv(diversity1, "Fungus Diversity non-transformed.csv")
write.csv(diversity2, "Fungus Diversity rarefied.csv")

plot_richness(physeq, measures=c("Shannon", "Simpson"))
plot_richness(ps.rarefied, measures=c("Shannon", "Simpson"))

df<-ordinate(physeq, method="PCoA", distance="bray", justDF=T)
write.csv(df$vectors, "Fungus  sample PCoA.csv")
ordination = ordinate(physeq, method="PCoA", distance="bray")
sample_data(physeq)['sample_id'] <- row.names(sample_data(physeq)) 
plot_ordination(physeq, ordination, "samples", label="sample_id") + theme_minimal()

df<-ordinate(ps.rarefied, method="PCoA", distance="bray", justDF=T)
write.csv(df$vectors, "Fungus rarefied sample PCoA.csv")
ordination = ordinate(ps.rarefied, method="PCoA", distance="bray")
sample_data(ps.rarefied)['sample_id'] <- row.names(sample_data(ps.rarefied)) 
plot_ordination(ps.rarefied, ordination, "samples", label="sample_id") + theme_minimal()
########################



#BACTERIA###############
data<-read.csv("Taxa level Correlations UPDATED/BACTERIA ASVs.csv", row.names = "ASV_ID") #Bacteria reps averaged
otumat<-data[1:48]
otumat[] <-lapply(otumat, as.integer)
otumat<-as.matrix(otumat)


taxmat<-data[49:54]
taxmat<-as.matrix(taxmat)
#######################

#make into phyloseq object
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)
physeq


#Rarefaction
sample_sums(physeq)
#rarecurve(t(otumat), step=50, cex=0.5)
ps.rarefied = rarefy_even_depth(physeq, rngseed=1, sample.size=0.9*min(sample_sums(physeq)), replace=F)
write.csv(otu_table(ps.rarefied), "Bacteria rarefied counts.csv")
plot_bar(physeq, fill = "phylum")
plot_bar(ps.rarefied, fill = "phylum")
otu_matrix <- as(otu_table(physeq), "matrix")
rarecurve(t(otu_matrix), step = 100, label = T)

diversity1<-estimate_richness(physeq)
diversity2<-estimate_richness(ps.rarefied)
write.csv(diversity1, "Bacteria Diversity non-transformed.csv")
write.csv(diversity2, "Bacteria Diversity rarefied.csv")

plot_richness(physeq, measures=c("Shannon", "Simpson"))
plot_richness(ps.rarefied, measures=c("Shannon", "Simpson"))

df<-ordinate(physeq, method="PCoA", distance="bray", justDF=T)
write.csv(df$vectors, "Bacteria  sample PCoA.csv")
ordination = ordinate(physeq, method="PCoA", distance="bray")
sample_data(physeq)['sample_id'] <- row.names(sample_data(physeq)) 
plot_ordination(physeq, ordination, "samples", label="sample_id") + theme_minimal()

df<-ordinate(ps.rarefied, method="PCoA", distance="bray", justDF=T)
write.csv(df$vectors, "Bacteria rarefied sample PCoA.csv")
ordination = ordinate(ps.rarefied, method="PCoA", distance="bray")
sample_data(ps.rarefied)['sample_id'] <- row.names(sample_data(ps.rarefied)) 
plot_ordination(ps.rarefied, ordination, "samples", label="sample_id") + theme_minimal()
