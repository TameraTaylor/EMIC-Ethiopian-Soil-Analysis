setwd("C:/Users/small/Documents/Brady lab/Soil Paper")

library(reshape2)
library(pheatmap)
library(RColorBrewer)


data <- read.csv("FullSoilV3.2.csv")





########10-3-2024 Without Striga info

data2<-data[-c(30),] #remove outlier E30
row.names(data2)<-data2[,1]
colnames(data2)
data3<-data2[, c(6:37)] #physio-chem
data3m<-data.matrix(data3)

#Row annotation of soil texture added May 2026
row_ann<-as.data.frame(data2[,c(1,38)])
row.names(row_ann)<-data2$ID
row_ann<-row_ann[-c(1)]
row_ann$USDA_Texture<-as.factor(row_ann$USDA_Texture)
palette_ramp <- colorRampPalette(c("white", "black"))(9)
print(palette_ramp)
ann_colors = scale_colour_gradient(low = "white", high = "black", breaks=9)
ann_colors = list(USDA_Texture=c("CLAY"="#FFFFFF", "CLAY LOAM"="#DFDFDF", "LOAM"= "#BFBFBF",  "SANDY CLAY"= "#9F9F9F","SANDY CLAY LOAM"= "#7F7F7F", "SILT"= "#5F5F5F",  "SILT LOAM"= "#3F3F3F", "SILTY CLAY"= "#1F1F1F", "SILTY CLAY LOAM"=  "#000000"))

col_ann<-read.csv("Physicochem annotations.csv")
col_ann$Type<-as.factor(col_ann$Type)
col_ann<-col_ann[-c(33:36),]
rownames(col_ann)<-col_ann[,1]
col_ann<-col_ann[-c(1)]
gp_row = which(diff(as.numeric(factor(col_ann$Type)))!=0)

pheatmap(data3m, cluster_cols = FALSE, cluster_rows = TRUE, clustering_distance_rows = "euclidean",
            scale="column", cutree_rows=6, border_color = "gray80",
            gaps_col = gp_row, annotation_row = row_ann, annotation_colors = ann_colors)
clusters <- cutree(p$tree_row, k = 6)
data_with_clusters <- cbind(cluster = clusters, data3m)
data_with_clusters<-as.data.frame(data_with_clusters)
cluster_means <- aggregate(. ~ cluster, data = data_with_clusters, FUN = mean)
write.csv(cluster_means, "physiochem cluster means.csv")

## ONLY STRIGA DATA
data2<-data[-c(30, 15, 48:50),] 
row.names(data2)<-data2[,1]
colnames(data2)
data3<-data2[, c(43:46)] #striga
data3m<-data.matrix(data3)
