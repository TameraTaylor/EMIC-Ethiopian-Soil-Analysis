library(tidyr)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(openxlsx)
library(ggplot2)
library(ggrepel)

setwd("C:/Users/small/Documents/Brady lab/Soil Paper/Luisa ASVs/Taxa level Correlations UPDATED/Imputed CLR Unknowns Grouped")

CLRFun<-read.csv("CLR_Fungus.csv", row.names = "X")
CLRBact<-read.csv("CLR_Bacteria.csv", row.names = "X")
FunTax<-read.csv("../FUNGUS_ASVS_TaxFilled.csv")[,c(1, 50:55)]
BactTax<-read.csv("../BACTERIA_ASVS_TaxFilled.csv")[,c(1, 50:55)]


########################PCA
#use imputed & clr transformed ASV counts to run analysis, then compare to taxanomic info of most informative things
pcaAnalysis<-function(df, tax, filename){
  pca_result <- prcomp(t(df), scale. = FALSE, center = TRUE) #prcomp cant handle NAs
  pc_scores <- pca_result$x
  loadings <- as.data.frame(pca_result$rotation)
  loadings$ASV_ID <- rownames(loadings)
  loadings_with_tax <- left_join(loadings, tax, by = "ASV_ID")
  row.names(loadings_with_tax)<-loadings_with_tax$ASV_ID
  top_PC1 <- loadings_with_tax %>% 
    arrange(desc(abs(PC1))) %>% 
    head(10)
  summary<-summary(pca_result)
  importance<-summary$importance
  
  wb <- createWorkbook()
  addWorksheet(wb, "scores")
  addWorksheet(wb, "loadings")
  addWorksheet(wb, "importance")
  addWorksheet(wb, "top10")
  writeData(wb, "scores", pc_scores, rowNames=TRUE)
  writeData(wb, "loadings", loadings_with_tax, rowNames=TRUE)
  writeData(wb, "importance", importance, rowNames=TRUE)
  writeData(wb, "top10", top_PC1, rowNames=TRUE)
  saveWorkbook(wb, paste(filename, ".xlsx", sep=""), overwrite = TRUE)
  
  pdf(file= paste(filename, ".pdf", sep=""), width=4, height=4)
  p<-autoplot(pca_result, data = t(df), loadings = FALSE) + theme_bw() + 
    geom_text_repel(aes(label= rownames(pc_scores)), max.overlaps = Inf,  alpha = 0.7)
  print(p)
  dev.off()
  
}


pcaAnalysis(CLRFun, FunTax, "PCA_Fungus_scaleFALSE")
pcaAnalysis(CLRBact, BactTax, "PCA_Bacteria_scaleFALSE")







########################################PCoA
# PCoA Analysis using imputed & clr transformed ASV counts
# PCoA is based on distance matrices (Euclidean distance for clr-transformed data)

pcoaAnalysis <- function(df, tax, filename){
  # Calculate Euclidean distance matrix (appropriate for clr-transformed data)
  dist_matrix <- dist(t(df), method = "euclidean")
  pcoa_result <- cmdscale(dist_matrix, k = ncol(df) - 1, eig = TRUE)
  
  # Extract scores (principal coordinates)
  pc_scores <- pcoa_result$points
  colnames(pc_scores) <- paste0("PCo", 1:ncol(pc_scores))
  
  # Calculate eigenvalues and variance explained
  eigenvalues <- pcoa_result$eig
  variance_explained <- eigenvalues / sum(abs(eigenvalues)) * 100
  
  # PCoA doesn't have direct loadings like PCA. Instead calculate correlations between original variables and PCo axes
  loadings <- cor(t(df), pc_scores, use = "complete.obs")
  colnames(loadings) <- paste0("PCo", 1:ncol(pc_scores))
  loadings <- as.data.frame(loadings)
  loadings$ASV <- rownames(loadings)
  tax$ASV <- rownames(tax)
  loadings_with_tax <- left_join(loadings, tax, by = "ASV")
  rownames(loadings_with_tax) <- loadings_with_tax$ASV
  
  # Save results to Excel
  wb <- createWorkbook()
  addWorksheet(wb, "scores")
  addWorksheet(wb, "loadings")
  addWorksheet(wb, "variance")
  writeData(wb, "scores", pc_scores, rowNames = TRUE)
  writeData(wb, "loadings", loadings_with_tax, rowNames = TRUE)
  writeData(wb, "variance", data.frame(
    Axis = paste0("PCo", 1:length(variance_explained)),
    Eigenvalue = eigenvalues,
    Variance_Percent = variance_explained,
    Cumulative_Variance = cumsum(variance_explained)
  ))
  saveWorkbook(wb, paste(filename, ".xlsx", sep = ""), overwrite = TRUE)
  
  # Create PCoA plot

  plot_data <- as.data.frame(pc_scores[, 1:2])
  colnames(plot_data) <- c("PCo1", "PCo2")
  p <- ggplot(plot_data, aes(x = PCo1, y = PCo2)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_text_repel(aes(label = rownames(plot_data)), size = 3) +
    labs(
      title = filename,
      x = paste0("PCo1 (", round(variance_explained[1], 2), "%)"),
      y = paste0("PCo2 (", round(variance_explained[2], 2), "%)")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )

  dev.new() 
  print(p)


  
  # Save summary statistics
  sink(paste(filename, ".txt", sep = ""))
  cat("PCoA Analysis Results\n")
  cat("=====================\n\n")
  cat("Eigenvalues and Variance Explained:\n")
  print(data.frame(
    Axis = paste0("PCo", 1:length(variance_explained)),
    Eigenvalue = eigenvalues,
    Variance_Percent = variance_explained,
    Cumulative_Variance = cumsum(variance_explained)
  ))
  cat("\n\nGower's distance (from eigenvalues):\n")
  print(summary(pcoa_result))
  sink()
}



pcoaAnalysis(CLRFun, FunTax, "PCoA_Fungus")
pcoaAnalysis(CLRBact, BactTax, "PCoA_Bacteria")

# Identify top ASVs contributing to PCo1 (highest correlation with axis)
# For Fungi
dist_matrix_fun <- dist(t(CLRFun), method = "euclidean")
pcoa_res_fun <- cmdscale(dist_matrix_fun, k = ncol(CLRFun) - 1, eig = TRUE)
loadings_fun <- cor(t(CLRFun), pcoa_res_fun$points, use = "complete.obs")
colnames(loadings_fun) <- paste0("PCo", 1:ncol(pcoa_res_fun$points))
loadings_fun_df <- as.data.frame(loadings_fun)
loadings_fun_df$ASV <- rownames(loadings_fun_df)
loadings_fun_df <- left_join(loadings_fun_df, FunTax, by = "ASV")
rownames(loadings_fun_df) <- loadings_fun_df$ASV

top_PC1_fungi <- loadings_fun_df %>% 
  arrange(desc(abs(PCo1))) %>% 
  head(10)
View(top_PC1_fungi)

# For Bacteria
dist_matrix_bact <- dist(t(CLRBact), method = "euclidean")
pcoa_res_bact <- cmdscale(dist_matrix_bact, k = ncol(CLRBact) - 1, eig = TRUE)
loadings_bact <- cor(t(CLRBact), pcoa_res_bact$points, use = "complete.obs")
colnames(loadings_bact) <- paste0("PCo", 1:ncol(pcoa_res_bact$points))
loadings_bact_df <- as.data.frame(loadings_bact)
loadings_bact_df$ASV <- rownames(loadings_bact_df)
loadings_bact_df <- left_join(loadings_bact_df, BactTax, by = "ASV")
rownames(loadings_bact_df) <- loadings_bact_df$ASV

top_PC1_bacteria <- loadings_bact_df %>% 
  arrange(desc(abs(PCo1))) %>% 
  head(20)
View(top_PC1_bacteria)