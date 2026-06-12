#Functions to live outside of the reactive codes

#ETHIOPIA BASE MAP

# MapBase<-  ggplot(data=r4) +
#     geom_raster(aes(x = x, y = y, fill = KG3_CRUTS32_Hist_8110raster)) +
#     geom_sf(data=ethiopia_regions_sf, fill=NA, linewidth=.8)+
#     scale_fill_manual(values=rev(mypalette))+ 
#     theme_void() 


#ETHIOPIA INTERACTIVE MAP
plotlymap <- function(pointColor, name){
   data1 <- data %>%
    mutate(txt=paste("ID: ", ID, "\nlongitude: ", lon, " ", "\nlatitude: ", lat, "\nSite: ", site, "\n", pointColor))
  unit <- SoilData_units[,2][SoilData_units$Name == name]
  if (!is.na (unit) & unit!="") {name <- paste(name, unit)}
  ggplotly(MapBase +
             geom_point(data = data1, aes(x = lon, y = lat, color = pointColor, text=txt), size = 3, alpha = 0.3) +
             scale_color_viridis(option="F", name = name) + labs(color=name)+
             coord_map()  , tooltip = "text", source = "Bubble Map")
}

#ETHIOPIA INTERACTIVE NO-SCALE MAP
noscalemap <- function(){
  ggplotly(MapBase +
             geom_point(data = data, aes(x = lon, y = lat, text=ID), color = "black", size = 3, alpha = 0.3) +
             coord_map(),  tooltip = "text", source = "No Scale Map")
}

################################
 #COMPOSITIONS
##############################
physiopie<-function(){
  id <- if (length(chosenid()) > 0) chosenid() else defaultID
    soilcomp1 <- c(data$organic, data$Carbonated_lime, data$Clay, data$Silt, data$Sand)[data$ID == id]
    soilcomp1[is.na(soilcomp1)] = 0
  pie(soilcomp1, labels = c("Organic matter", "Carbonated lime", "Clay (<2 um)", "Silt (2-50 um)", "Sand (>50 um)"), border="white", col = rainbow(length(soilcomp1)),
      main = "Physical Attributes")
}

nut1fun<- function(){
  mgkgchoices<-c("ID", "N", "S_avail",	"S_total",	"P_avail",	"P_total",	"K_avail",	"Mg_avail", "Na_avail")
  mgkg<-data[mgkgchoices]
  id <- if (length(chosenid()) > 0) chosenid() else defaultID
    mgkg<-mgkg[mgkg$ID == id,]#comma after logical indicates rows
    mgkg[is.na(mgkg)] = 0
  Meltedmgkg<-reshape2::melt(mgkg, id=c("ID"), variable.name="Variable", value.name="Value")
  ggplot(Meltedmgkg, aes(x=Variable,y=Value, fill=Variable)) +geom_bar(stat = 'identity') + labs(y="mg/kg", x="") + theme(legend.position='none') + 
    ggtitle("Nutrient Composition") 
}

nut2fun<- function(){
  mmkgchoices<-c("ID", "K_total",	"Ca_total",	"Mg_total",	"Na_total")
  mmkg<-data[mmkgchoices]
  id <- if (length(chosenid()) > 0) chosenid() else defaultID
    mmkg<-mmkg[mmkg$ID == id,]#comma after logical indicates rows
    mmkg[is.na(mmkg)] = 0
  Meltedmmkg<-reshape2::melt(mmkg, id=c("ID"), variable.name="Variable", value.name="Value")
  ggplot(Meltedmmkg, aes(x=Variable,y=Value, fill=Variable)) +geom_bar(stat = 'identity') + labs(y="mmol/kg", x="") + theme(legend.position='none')
}

nut3fun<-function(){
  mmlchoices<-c("ID", "Ca_avail")
  mml<-data[mmlchoices]
  id <- if (length(chosenid()) > 0) chosenid() else defaultID
    mml<-mml[mml$ID == id,]#comma after logical indicates rows
    mml[is.na(mml)] = 0
  Meltedmml<-reshape2::melt(mml, id=c("ID"), variable.name="Variable", value.name="Value")
  ggplot(Meltedmml, aes(x=Variable,y=Value, fill=Variable)) +geom_bar(stat = 'identity') + labs(y="mmol/l", x="") + theme(legend.position='none')
}

#microbe pie
microbepie<-function(df, titlename, idchoice){
  id <- if (length(chosenid()) > 0) chosenid() else defaultID
  df$phylum[is.na(df$phylum)] <- "Unknown"
  microbes <- df %>%
    group_by(phylum) %>%
    summarise(across(starts_with("E"), sum, na.rm = TRUE))

  phylum_counts <- setNames(microbes[[id]], microbes$phylum)
  phylum_counts <- sort(phylum_counts, decreasing = TRUE)

  top10 <- head(phylum_counts, 10)
  other <- sum(phylum_counts[-(1:10)])
  pie_data <- data.frame(
    phylum = as.factor(c(names(top10), "Other")),
    abundance = c(top10, other)
  )
  pie_data$phylum <- factor(pie_data$phylum, levels = pie_data$phylum)
  pie_data$Percentage <- round(pie_data$abundance / sum(pie_data$abundance) * 100, 1)

  
  ggplot(pie_data, aes(x = "", y = abundance, fill = phylum)) +
    geom_bar(stat = "identity", width = 1, color="white") +
    scale_fill_manual(values=mycolors)+
    coord_polar("y") + 
    theme_void()+
    # geom_text(aes(label=paste(Percentage, "%")),
    #           position=position_stack(vjust=0.5),
    #           size=5)+
    geom_label_repel(data = pie_data,
                     aes(label = paste0(Percentage, "%")),
                     size = 4.5, nudge_x = 1, show.legend = FALSE)+
    ggtitle(titlename)
    
}


#############################
  #CORRELATIONS
###########################

#CORR
CorrelationFun<-function(CorrType){
  FilteredDF<-rdf()
  UnmeltedFDF <- dcast(FilteredDF, ID + lon + lat ~ DataType, value.var="Value")
  xx<-UnmeltedFDF[,-(1:3)] #remove ID, lon & lat columns
  xx[is.na(xx)] = 0 #change NAs to 0s. Need to think of a more appropriate solution for this as 49 & 50 are just missing microbe data ad aren't 0s.
  xx<-remove_constant(xx, na.rm = FALSE, quiet = TRUE)
  C <- cor(xx,  use="pairwise.complete.obs", method=paste(CorrType))
  cortest<-corr.test(xx,method = paste(CorrType))$p 
  corrplot(C,type='lower',  
           method = 'color', 
           diag = FALSE, 
           p.mat = cortest, 
           insig = 'label_sig', 
           sig.level = c(0.001, 0.01, 0.05), 
           order="hclust", 
           tl.col="black")
} 
# CorTableFun<-function(CorrType){
#     FilteredDF<-rdf()
#     UnmeltedFDF <- dcast(FilteredDF, ID + lon + lat ~ DataType, value.var="Value")
#     xx<-UnmeltedFDF[,-(1:3)]
#     C <- cor(xx,  use="pairwise.complete.obs", method=paste(CorrType))
#     C
#   }

#PCA
PCAFun<-function(){
  FilteredDF<-rdf()
  UnmeltedFDF <- dcast(FilteredDF, ID + lon + lat ~ DataType, value.var="Value")
  row.names(UnmeltedFDF)<-UnmeltedFDF[,1]
  DataMatPRC<-data.matrix(UnmeltedFDF[,-(1:3)]) #make the values except ID, lon & lat into a matrix
  DataMatPRC[is.na(DataMatPRC)] = 0 #change NAs to 0s. Need to think of a more appropriate solution for this as 49 & 50 are just missing microbe data ad aren't 0s.
  DataMatPRC<-remove_constant(DataMatPRC, na.rm = FALSE, quiet = TRUE)
  soilpca<-prcomp(DataMatPRC, scale = TRUE)
  autoplot(soilpca, data=UnmeltedFDF, label=TRUE, shape=FALSE)
}

#HEATMAP
HeatmapFun<-function(){
  FilteredDF<-rdf()
  UnmeltedFDF <- dcast(FilteredDF, ID + lon + lat ~ DataType, value.var="Value")
  row.names(UnmeltedFDF)<-UnmeltedFDF[,1]
  UnmeltedFDF<-UnmeltedFDF[,-(2:3)] #remove lon & lat
  DataMat<-data.matrix(UnmeltedFDF[,-(1)]) #do not include ID in matrix
  DataMat<-remove_constant(DataMat, na.rm = FALSE, quiet = TRUE)
  print(head(DataMat))
  return(DataMat)
}
#   p<-pheatmap(DataMat, cluster_cols = TRUE, cluster_rows = TRUE, clustering_distance_rows = "manhattan",
#               scale="column", border_color = "gray80", cellwidth = 10, cellheight = 10)
#   dev.off()
#   p
# }

#PAIRS CORRELATION
panel.rr <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  rr <- summary(lm(x~y))$r.squared # Remove abs function if desired
  txt <- paste0(prefix, format(c(rr, 0.123456789), digits = digits)[1])
  text(0.5, 0.5, txt) # Resize the text by level of correlation
}
PairsFun<-function(){
  FilteredDF<-rdf()
  UnmeltedFDF <- dcast(FilteredDF, ID + lon + lat ~ DataType, value.var="Value")
  xx<-UnmeltedFDF[,-(1:3)] #remove ID, lon & lat columns
  xx[is.na(xx)] = 0
  pairs(xx, lower.panel=reg, upper.panel=panel.rr) #graphics package
}
reg <- function(x, y, ...) {
  points(x,y, ...)
  abline(lm(y~x)) 
}  #makes regression line


#############################
  #MODELING
############################  

DecisionTree<- function(){
  tree <- rpart(paste(input$ModDVariable, "~ ."), data=Modrdf(), control=rpart.control(cp=.0000000001))
  prp(tree,
      faclen=0, #use full names for factor labels
      extra=1, #display number of obs. for each terminal node
      roundint=F, #don't round to integers in output
      digits=5) #display 5 decimal places in output
}

treeimportanceFUN<-function(variable){
  set.seed(71)
  modd<-variable
  formula<-paste(modd, "~ .")
  formula<-as.formula(formula)
  rf <-randomForest(formula, data=Modrdf(), importance=TRUE, ntree=500)  #Evaluate variable importance 
  #importance(rf)
  #varImpPlot(rf)
  imp = as.data.frame(importance(rf))
  imp = cbind(vars=rownames(imp), imp)
  imp = imp[order(imp$"%IncMSE"),]
  imp$vars = factor(imp$vars, levels=unique(imp$vars))
  imp %>% 
    pivot_longer(cols=matches("Inc")) %>% 
    ggplot(aes(value, vars)) +
    geom_col() +
    geom_text(aes(label=round(value), x=0.5*value), size=3, colour="white") +
    facet_grid(. ~ name, scales="free_x") +
    scale_x_continuous(expand=expansion(c(0,0.04))) +
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          axis.title=element_blank())
}



#############
#MICROBE
#############

MicrobeAggregate<- function(Taxa1, Taxa2){
  if(Taxa1 == "Bacteria"){
    Microbe_agg <- BacteriaData}
  if(Taxa1 == "Fungi"){
    Microbe_agg<- FungiData}
  Microbe_agg[is.na(Microbe_agg)] <- "Unknown"
  Microbe_agg<- dplyr::select(Microbe_agg, all_of(input$ID2),input$Taxa2)
  if(Taxa2 != "ASV"){
  Microbe_agg<- Microbe_agg %>% group_by(!!rlang::sym(Taxa2)) %>%
    summarise(across(starts_with("E"), sum, na.rm = TRUE)) %>% #sum across Soil Samples
    as.data.frame() }
  else{
  Microbe_agg<-dplyr::as_tibble(Microbe_agg, rownames = "ASV") %>% as.data.frame()
  }
  Microbe_agg
}

#microbe abundance barplot for all samples
MicrobeBarPlot<-function(Taxa1, Taxa2){
  df<-MicrobeAggregate(Taxa1, Taxa2)
  row.names(df)<-df[,1]

  total_sums <- colSums(df %>% dplyr::select(starts_with("E")), na.rm = TRUE)
  proportions <- df %>%
    mutate(across(starts_with("E"), ~ . / total_sums[cur_column()])) #%>% dplyr::select(1:49)
  total_taxa <- rowSums(df %>% dplyr::select(starts_with("E")), na.rm = TRUE)
  top10 <- names(sort(total_taxa, decreasing = TRUE)[1:10])
  top10_data <- proportions %>% filter(rownames(.) %in% top10)
  other_data <- proportions %>% filter(!rownames(.) %in% top10)
  other_row <- colSums(other_data %>% dplyr::select(starts_with("E")), na.rm = TRUE)
  plot_data <- rbind(top10_data[-1], Other = other_row)
  plot_data <- plot_data %>%
    tibble::rownames_to_column(var = "Taxon")
  df_bar_long <- plot_data %>%
    pivot_longer(cols = -Taxon, names_to = "Sample", values_to = "Abundance")
  df_bar_long$Taxon <- factor(df_bar_long$Taxon, levels = c(setdiff(unique(df_bar_long$Taxon), "Other"), "Other"))
  
  ggplot(df_bar_long, aes(x = Sample, y = Abundance, fill = Taxon)) +
    geom_bar(position = "fill", stat = "identity") +
    labs(x = "Soil ID", y = "Percentage") + scale_fill_manual(values=mycolors)+
    theme_bw() + scale_y_continuous(labels = scales::percent_format()) +
    theme(legend.position = "top", axis.text.x=element_text(angle=60,hjust=1)) + labs(fill=Taxa2) 
}

MicrobeCLR<-function(Taxa1, Taxa2){
  df<-MicrobeAggregate(Taxa1, Taxa2)
  row.names(df)<-df[,1]
  set.seed(71)
  CLR_df<-clr_lite(df[,-1], samples_are = "cols", method="logunif", replicates=1000)
  CLR_df
}


create_microbe_df <- function(clr_df, traits) {
  print(traits)
  df1 <- t(clr_df) %>% as.data.frame()
  df1$ID <- rownames(df1)
  df2 <- dplyr::select(data, ID, all_of(traits))
  merge(df1, df2, by = "ID")
}


microbe_corr <- function(dataT, method) {
  tax_cols <- 2:(ncol(dataT) - length(input$variable2))
  trait_cols <- (ncol(dataT) - length(input$variable2) + 1):ncol(dataT)
  dataTax<-as.matrix(dataT[, tax_cols])
  dataTrait<-as.matrix(dataT[, trait_cols])
  corr.test(dataTax, dataTrait, 
            method = method, adjust = "bonferroni")
}

microbeHeatmap<-function(){
  cor_vals <- cor_results()$r %>% base::as.matrix()
  significance<-cor_results()$p %>% base::as.matrix()
  significance <- cut(significance, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                      labels = c("***", "**", "*", ""), right = FALSE)
  print(head(cor_vals))
  print(head(significance))
  rg<-max(abs(cor_vals))
  p<-pheatmap(cor_vals, display_numbers=significance, main = "Microbe-Trait Correlation", fontsize = 10, fontsize_number = 10, color=colorRampPalette(c("red", "white", "blue"))(30), breaks=seq(-rg, rg, length.out=30) )
  dev.off()
  p
}


