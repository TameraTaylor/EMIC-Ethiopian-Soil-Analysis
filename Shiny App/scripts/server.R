server <- function(input, output, session) {
  #give access to functions
  source("scripts/Functions.R", local=TRUE)

  
  #intro photo
  output$photo1 <- renderImage({
    list(src = "www/field_photo.png", height=400
    )}, deleteFile = FALSE)

  
  ###########################
  #MAPS
  ############################
  
  # Main Map
  output$mainmap <- renderPlotly({
    c <- data[[input$colorDet]]
      plotlymap( c, input$colorDet)
  })

  
  ######################################
               #COMPOSITIONS 
  #####################################
  
  # For the default pie charts
  defaultID <- "E01"  
  
  #Display the name of the sample that is picked  
  output$e <- renderText({
    d <- event_data("plotly_click", source = "No Scale Map")
    id <- data$ID[data$lon == d$x & data$lat == d$y]
    if (!is.null(d) & length(id) > 0){
      id
    }else{
      defaultID
    }
  }) 
  
  chosenid<-reactive({d <- event_data("plotly_click", source = "No Scale Map")
  id <- data$ID[data$lon == d$x & data$lat == d$y]
  id
  })
  
  #Map for nutrient selection
  output$piechartmap <- renderPlotly({
    noscalemap()
  })
  
  output$pie1 <- renderPlot({
    physiopie()
  })

  output$nutrients1 <- renderPlot({
    nut1fun()
  },height=300)
  
  output$nutrients2 <- renderPlot({
    nut2fun()
  },height=300)
  
  output$nutrients3 <- renderPlot({
    nut3fun()
  },height=300)
  
  #Microbes
  output$pie3 <- renderPlot({
    microbepie(BacteriaData, "Bacterial Phyla", idchoice)
  })
  output$pie4<-renderPlot({
    microbepie(FungiData, "Fungal Phyla", idchoice)
  })
  
  #Downloading
  output$downloadpie1 <- downloadHandler(
    filename = function() {
      if (length(chosenid()) > 0){
        paste(chosenid(), "-physicalcomp.pdf", sep="")
      }
      else{
        paste(defaultID, "-physicalcomp.pdf", sep="")
      }
    },
    content = function(file) {
      pdf(file)
      print(physiopie())
      dev.off()
    })
  
  output$downloadpie2 <- downloadHandler(
    filename = function() {
      if (length(chosenid()) > 0){
        paste(chosenid(), "-nutrients.pdf", sep="")
      }
      else{
        paste(defaultID, "-nutrients.pdf", sep="")
      }
    },
    content = function(file) {
      pdf(file)
      arrangeGrob( 
        print(nut1fun()), 
        print(nut2fun()),
        print(nut3fun()), ncol = 3)
      dev.off()
    })
  
  output$downloadpie3 <- downloadHandler(
    filename = function() {
      if (length(chosenid()) > 0){
        paste(chosenid(), "-microbes.pdf", sep="")
      }
      else{
        paste(defaultID, "-microbes.pdf", sep="")
      }
    },
    content = function(file) {
      pdf(file)
      print(microbepie(BacteriaData, "Bacterial Phyla"))  
      print(microbepie(FungiData, "Fungal Phyla"))
      dev.off()
    })
  

  
  
  ######################################################
                 #CORRELATIONS 
  #####################################################

  output$CorrTypeSelect <- renderUI({ #only displays this option when correlation is chosen as type of analysis
    if (input$graph == "Correlation Matrix") {
    selectInput(inputId="CorrType", "Choose correlation type:",
                choices=c("pearson", "spearman"), multiple=FALSE)
    }
  })

  mdf<-reactive({
    Melteddata<-reshape2::melt(data, id=c("ID", "lon", "lat"), variable.name="DataType", value.name="Value")
    return(Melteddata)
  })
  
  rdf<-eventReactive(input$corrGO,{
    if (!("All" %in% input$ID)) {
    FilteredDF<- mdf()%>%
      filter(ID %in% input$ID)}
    else {FilteredDF<- mdf()}
    FilteredDF<-FilteredDF%>% 
        filter(DataType %in% input$variable)
    FilteredDF$Value<-as.numeric(FilteredDF$Value)
    return(FilteredDF)
  })
  
 eventReactive(input$corrGO,{
   updateSelectInput(session=session,inputId="ID",label = "Choose the samples you want to see:")
   updateSelectInput(inputId = "variable")
   updateSelectInput(session=session, inputId="graph", "Choose a type of Analysis", choices = c("Correlation Matrix",
                                                                                                "Pairs Scatter Plots",
                                                                                                "Principal Component Analysis",
                                                                                                "Heatmap"))
  })  
 
  output$plots <- renderPlot({input$corrGo
    graph<-isolate(input$graph)
    if(graph == "Correlation Matrix") {
     CorrelationFun(input$CorrType)
    }
    else if (graph == "Pairs Scatter Plots") {
     PairsFun()
    }
    else if (graph == "Principal Component Analysis") {
      PCAFun()
    } #cannot handle NAs from microbe data
    else if (graph == "Heatmap") {
     DataMat<-HeatmapFun()
      p<-pheatmap(DataMat, cluster_cols = TRUE, cluster_rows = TRUE, clustering_distance_rows = "manhattan",
                  scale="column", border_color = "gray80", cellwidth = 10, cellheight = 10)
      dev.off()
      p
    }
  })
  
 # output$CorText<-renderTable({
 #    CorTableFun(input$CorrType)
 #     })
 
  output$stattext <- renderText({
    req(input$corrGO > 0)
    isolate({
    if (input$graph == "Correlation Matrix") {
      "Correlelogram displaying the correlation values between all variable chosen as colors. ***, **, * represent p-values less than or equal to 0.001, 0.01, and 0.05 respectively."
    }
    else if (input$graph == "Pairs Scatter Plots"){
      "Scatterplots for each pair of variables displayed with a linear regression line. The r squared value is displayed in the upper panels."
    }
    else if (input$graph == "Principal Component Analysis"){
      "Principal Component Analysis (PCA) is a way to reduce the data into two primary dimensions, assuming linear relationships, that are able to explain the most amount of the variation in the data. 
    It is useful for vizualizing how the samples relate to one another. Samples that cluster together are more similar to each other than those that are far apart."
    }
    else if (input$graph == "Heatmap"){
      "This heatmap is displaying the scaled values and has ordered the selected variable and samples by hierarchical clustering. 
    Items closer to eachother in the tree behave more similarly."
    }
    })
  })


  output$downloadcorrtab <- downloadHandler(
    filename = function() {
      if (input$graph == "Correlation Matrix"){
        paste("correlation_", format(strptime(Sys.time(), "%Y-%m-%d %H:%M:%S"), "%Y-%m-%d-%H-%M-%S"), ".pdf", sep="")
      }
      else if (input$graph == "Pairs Scatter Plots"){
        paste("scatter_", format(strptime(Sys.time(), "%Y-%m-%d %H:%M:%S"), "%Y-%m-%d-%H-%M-%S"), ".pdf", sep="")
      }
      else if (input$graph == "Principal Component Analysis"){
        paste("PCA_", format(strptime(Sys.time(), "%Y-%m-%d %H:%M:%S"), "%Y-%m-%d-%H-%M-%S"), ".pdf", sep="")
      } 
      else if (input$graph == "Heatmap"){
        paste("heatmap_", format(strptime(Sys.time(), "%Y-%m-%d %H:%M:%S"), "%Y-%m-%d-%H-%M-%S"), ".pdf", sep="")
      }
    },
    content = function(file) {
      if (input$graph == "Correlation Matrix"){
       pdf(file)
        arrangeGrob( 
        print(CorrelationFun()), 
        print(CorTableFun()), ncol = 2)
        dev.off()
      }
      else if (input$graph == "Pairs Scatter Plots"){
        pdf(file)
          print(PairsFun())
        dev.off()
      } 
      else if (input$graph == "Principal Component Analysis"){
        pdf(file)
        print(PCAFun())
        dev.off() 
      }
      else if (input$graph == "Heatmap"){
        DataMat <- HeatmapFun()
        n_rows <- nrow(DataMat)
        n_cols <- ncol(DataMat)
        pdf_width <- (n_cols * 10 + 200) / 72 # 72 pts per inch
        pdf_height <- (n_rows * 10 + 200) / 72
        p<-pheatmap(DataMat, cluster_cols = TRUE, cluster_rows = TRUE, clustering_distance_rows = "manhattan",
                    scale="column", border_color = "gray80", cellwidth = 10, cellheight = 10)

        pdf(file, width=pdf_width, height=pdf_height)
        print(p)
        dev.off() 
      }
    })

  
##################################################  
  #MODELING
##################################################
  
  
  Modrdf<-eventReactive(input$modelingGO,{
    targets<-c(input$ModIVariable,input$ModDVariable)
    Melteddata2<-reshape2::melt(data, id="ID", variable.name="DataType", value.name="Value") #NEW
    if (!("ALL" %in% input$ModID)) {
      FilteredDF<- Melteddata2%>%
        filter(ID %in% input$ModID)}
    else {FilteredDF<- Melteddata2}
    FilteredDF<-FilteredDF%>% 
      filter(DataType %in% targets)
    FilteredDF$Value<-as.numeric(FilteredDF$Value)
    Modrdf <- dcast(FilteredDF, ID ~ DataType, value.var="Value")
    Modrdf<- Modrdf[-c(1)] #exclude ID label
    Modrdf[is.na(Modrdf)] = 0
    return(Modrdf)
  })
  
  #output$forestdf<-renderTable(Modrdf())
  
  output$histogram<-renderPlot({
    hist(Modrdf()[[input$ModDVariable]], xlab = input$ModDVariable, main = paste("Histogram of" , input$ModDVariable))
  })#make note that Error: "invalid number of 'breaks'" will appear if dependent variable is only NA values
  
  
  
  output$tree <-renderPlot({
    DecisionTree()})
  
  output$treestats<-renderPrint({
    tree <- rpart(paste(input$ModDVariable, "~ ."), data=Modrdf(), control=rpart.control(cp=.0000000001))
    printcp(tree)
  })
  
  output$randomforest<-renderPrint({
    set.seed(71)
    modd<-input$ModDVariable
    formula<-paste(modd, "~ .")
    formula<-as.formula(formula)
    rf <-randomForest(formula, data=Modrdf(), importance=TRUE, ntree=500)
    print(rf)
  })
  
  
  output$treeimportance<-renderPlot({
    treeimportanceFUN(input$ModDVariable)
  })
    

  
  output$treetext<-renderText({
    "Decision trees and random forest model fits a training subset of the data to the model and then tests its 
  prediction accuracy against the unused rest of the dataset. With our limited number of samples, it was opted to use all the data in training 
  the model, therefore these models are for illustrative purposes only and do not have measures of significance."
  })
  
  
  
##########################################
  #MICROBES
###########################################

  
  aggregated_data <- eventReactive(input$MicrobeButton1, {
    MicrobeAggregate(input$Taxa1, input$Taxa2)
  })
  clr_data <- eventReactive(input$MicrobeButton1, {
    MicrobeCLR(input$Taxa1, input$Taxa2)
  })
  plot_bar <- eventReactive(input$MicrobeButton1, {
    MicrobeBarPlot(input$Taxa1, input$Taxa2)
  })
  
  output$TaxaAbundance <- renderPlot({
    req(plot_bar())
    plot_bar()
  })
  output$MicrobeBarDownload <- downloadHandler(
    filename = function() {
      paste(input$Taxa1, "-", input$Taxa2, "-Barplot.pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      print(plot_bar())
      dev.off()
    }
  )
  output$RawChosenTaxa <- DT::renderDataTable({
    req(aggregated_data())
    head(aggregated_data())
  }, options = list(dom = 't'))
  output$RawChosenTaxaDownload <- downloadHandler(
    filename = function() {
      paste(input$Taxa1, "-", input$Taxa2, "-RawAbundance.csv", sep = "")
    },
    content = function(file) {
      write.csv(aggregated_data(), file)
    }
  )
  output$CLRChosenTaxa <- renderUI({
    if (input$Taxa2 == "ASV") {
      HTML("<p>CLR transformation on the ASV data is too computationally intensive here...</p>")
    } else {
      div(style = 'overflow-x: scroll', DT::dataTableOutput("CLRTable"))
    }
  })
  output$CLRTable <- DT::renderDataTable({
    req(input$Taxa2 != "ASV")
    req(clr_data())
    head(clr_data())
  }, options = list(dom = 't'))
  output$CLRChosenTaxaDownload <- renderUI({
    if (input$Taxa2 != "ASV") {
      downloadButton("DownloadCLR", "Download Full Data")
    }
  })
  
  output$DownloadCLR <- downloadHandler(
    filename = function() {
      paste(input$Taxa1, "-", input$Taxa2, "-CLR.csv", sep = "")
    },
    content = function(file) {
      write.csv(clr_data(), file)
    }
  )
  
  
################
  #MICROBE ANALYSIS
##################
  merged_data <- eventReactive(input$MicrobeButton2, {
    req(input$variable2)
    create_microbe_df(clr_data(), input$variable2)
  })
  
  cor_results <- eventReactive(input$MicrobeButton2, {
    microbe_corr(merged_data(), input$CorrType2)
  })
  
  output$MicrobeChoices <- renderText({
    paste(
      "Summary of how you filtered your data:\n",
      "Samples chosen:", toString(input$ID2), "\n",
      "Taxa filtered at:", input$Taxa1, input$Taxa2
    )
  })
 
 # output$microbeAnalysis<-renderDataTable({
 #   req(input$variable2)
 #   merged_data()
 # })

 output$microbeCorr<-renderDataTable({
   cor_results()$r
 })
   
output$downloadcorrtab2<-downloadHandler(
  filename = function() {
      paste("microbe_correlation_", format(strptime(Sys.time(), "%Y-%m-%d %H:%M:%S"), "%Y-%m-%d-%H-%M-%S"), ".xlsx", sep="")
  },
  content = function(file) {
    wb <- createWorkbook()
    addWorksheet(wb, "P")
    addWorksheet(wb, "Padj")
    addWorksheet(wb, "R")
    writeData(wb, "P", cor_results()$p, rowNames=TRUE)
    writeData(wb, "Padj", cor_results()$p.adj, rowNames=TRUE)
    writeData(wb, "R", cor_results()$r, rowNames=TRUE)
    saveWorkbook(wb, file=file)},
  contentType = "file/xlsx"
)

output$MicrobeCorrelogram <- renderPlot({
 microbeHeatmap()
})

output$downloadMicroCorrPlot <- downloadHandler(
  filename = function() {
    paste("microbe_correlogram_", format(strptime(Sys.time(), "%Y-%m-%d %H:%M:%S"), "%Y-%m-%d-%H-%M-%S"), ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    print(microbeHeatmap())
    dev.off()
  }
)

}
