


setwd("C:/Users/small/Documents/Brady lab/Soil Paper")


library(ggplot2)


#need to center and scale data first
# Center and scale the dataset

data<-read.csv("FullSoilV7.csv")
str(data)
#log strigas are columns 42 and 44
physiochem <- as.matrix(data[,c(which(colnames(data)=="Mg_avail"):which(colnames(data)=="Clay_humus"))]) #removed agroclimate, need to recalc PCA
rownames(physiochem)<-data$ID
physiochem_df <- as.data.frame(physiochem)
physiochem_df<-physiochem_df[-c(30),]




#####################################

#11-9-2024
#using prcomp instead cuz used more often

#######################################
library(ggfortify)
summary(physiochem_df)
physiochem_df<-physiochem_df[,-19] #C.OS ratio had NAs, removed


pca_result <- prcomp(physiochem_df, scale. = TRUE, center = TRUE) #prcomp cant handle NAs
pc_scores <- pca_result$x
loadings <- pca_result$rotation
p<-autoplot(pca_result, data = physiochem_df,
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3, loadings.label.repel=T) + theme_bw()
d<-c(4, 5, 7,16,10, 11 ) #to only show some vectors
p$layers[[2]]$data[, d]
p$layers[[2]]$data[d]
p$layers[[2]]$data[d,]
p$layers[[2]]$data<-p$layers[[2]]$data[d,]
p$layers[[3]]$data<-p$layers[[3]]$data[d,]
p
write.csv(pc_scores, "C:/Users/small/Documents/Brady lab/Soil Paper/PCA for modeling/physiochem_pcscores.csv")
write.csv(loadings, "C:/Users/small/Documents/Brady lab/Soil Paper/PCA for modeling/physiochem_loadings.csv")
summary(pca_result)

