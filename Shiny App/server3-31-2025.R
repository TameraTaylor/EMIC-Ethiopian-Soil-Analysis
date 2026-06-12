server <- function(input, output, session) {
  #give access to functions
  source("./Functions.R", local=TRUE)

  
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
  },height=200, width=300)
  
  output$nutrients2 <- renderPlot({
    nut2fun()
  },height=200, width=300)
  
  output$nutrients3 <- renderPlot({
    nut3fun()
  },height=200, width=300)
  
  #Microbes
  output$pie3 <- renderPlot({
    microbepie()
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
      print(microbepie())
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
   updateSelectInput(session=session, inputId="graph", "Choose a type of Analysis", choices = c("Correlation Matrix",
                                                                                                "Pairs Scatter Plots",
                                                                                                "Principal Component Analysis",
                                                                                                "Heatmap"))
  })  
 
  output$plots <- renderPlot({
    if(input$graph == "Correlation Matrix") {
     CorrelationFun(input$CorrType)
    }
    else if (input$graph == "Pairs Scatter Plots") {
     PairsFun()
    }
    else if (input$graph == "Principal Component Analysis") {
      PCAFun()
    } #cannot handle NAs from microbe data
    else if (input$graph == "Heatmap") {
     HeatmapFun()
    }
  })
  
 output$CorText<-renderTable({
    CorTableFun(input$CorrType)
    })
  output$stattext <- renderPrint({
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


  output$downloadcorrtab <- downloadHandler(
    filename = function() {
      if (input$graph == "Correlation Matrix"){
        paste(format(strptime(Sys.time(), "%Y-%m-%d %H:%M:%S"), "%Y-%m-%d-%H-%M-%S"), "-correlation.pdf", sep="")
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
        pdf(file)
        print(HeatmapFun())
        dev.off() 
      }
    })

  
##################################################  
  #MODELING
##################################################
  
  # observe({
  #   updateSelectInput(session=session,inputId="ModID",label = "Choose the samples you want to model with:",
  #                     choices=c("ALL", data$ID))
  # }) #observe ID selection
  # observe({
  #   updateSelectInput(session=session,inputId="ModDVariable",label = "Choose the DEPENDENT variable you want to model with:",
  #                     choices=levels(mdf()$DataType))
  # }) #observe variable selection
  # observe({
  #   updateSelectInput(session=session,inputId="ModIVariable",label = "Choose the INDEPENDENT variables you want to model with:",
  #                     choices=levels(mdf()$DataType))
  # }) #observe variable selection
  
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
  

}
