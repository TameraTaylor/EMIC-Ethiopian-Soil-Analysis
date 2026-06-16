#separate ui for shiny app


sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Home", tabName = "home"),
    menuItem("Ethiopian Soil Summary",
             menuSubItem("Map View", tabName = "map"),
             menuSubItem("Individual Soil Compositions", tabName = "soilcomposition"),
             menuSubItem("Correlation Analysis", tabName = "correlationanalysis"),
             menuSubItem("Soil Modeling", tabName = "modeling")),
    menuItem("Ethiopian Soil Microbiome",
             menuSubItem("Data Browser & Filtering", tabName = "Data_Browser"), #just display the datatables of ASVs & taxonomy
             menuSubItem("Microbe Analysis", tabName="Microbe_Analysis")
             ) ) )




body <- dashboardBody(
  tabItems(
    tabItem(tabName = "home", 
            h2("Exploration of Ethiopian Soil Data (2010)"),
            h4("This data has been collected and analyzed in a collaborative effort between institutions as presented in 'Disentangling the importance of microbiological and physico-chemical properties of Ethiopian field soils for the Striga seed bank and sorghum infestation'.",
               "This application has been designed to offer means of exploratory analysis of these data to a wider audience.",
               "Use the side menu items to navigate through different aspects of the soil data and select options to view results in various ways.",
               "Future plans are to adapt the application to analyze user inputted data and provide treatment reccomendations for Striga infested soils."),
            imageOutput("photo1"), #also add logos of collaborators
            h5("Representative photo of sorghum field infested with Striga in Ethiopia. Photo courtesy: Dr. Taye Tessema"),
            h2("Sample Description"),
            h4("A total of 48 composite soil samples covering a trajectory of more than 1500 km were collected from naturally Striga infested sorghum growing agro ecological zones in Ethiopia in 2010.",
               "Soil samples were collected from naturally Striga infested sorghum fields in Amhara (Kemise, North Shewa, South and North Wollo Zones) and Tigray (West, Central and South zones) regions of Ethiopia.",
               "Sorghum fields with four categories (zero, low, medium and high) of Striga field infestation were randomly selected, and samples were taken from the top layer (0-20 cm) around the root zone of ",
               "a sorghum plant then combined together to form one composite sample per field.",
               "Physical, chemical, and microbiological compositions were analyzed for these samples."),
            h2("Aknowledgements"),
            h5("Tamera Taylor*1, Getahun Benti*2,3,4,†, Marcio F.A. Leite2,6, Luisa Arias-Giraldo2, Desalegn W. Etalo2,6, Sewunet Abera2,4, Lorenzo Lombard5, José G. Maciá-Vicente2,8, Stefan Sanow1, Dominika Rybka2, Thomas Mostert3, Einar Martinez de la Parte2, Daniel Yimer Legesse4, Urgesa Tsega Tulu4, Jiregna Daksa1,4,9, Roy van Doorn5, Raycenne R. Leite2, Taye Tessema4, Pedro W. Crous3, Dorota Kawa1,10,11, Eiko E. Kuramae2,3,12, Jos M. Raaijmakers2,3*, Siobhan M. Brady1,13*", br(),  
                "1 Department of Plant Biology and Genome Center, University of California, Davis, Davis, CA 95616", br(), 
                "2 Department of Microbial Ecology, Netherlands Institute of Ecology (NIOO -KNAW), Droevendaalsesteeg 10, 6708 PB Wageningen, The Netherlands", br(), 
                "3 Institute of Biology, Leiden University, 2333 BE Leiden, The Netherlands", br(), 
                "4 Ethiopian Institute of Agricultural Research, Holeta, Ethiopia", br(), 
                "5 Westerdijk Fungal Biodiversity Institute (WI-KNAW), Uppsalalaan 8, 3584 CT Utrecht, The Netherlands", br(), 
                "6 Plant Hormone Biology Group, Swammerdam Institute of Life Science, University of Amsterdam, Amsterdam, Netherlands", br(), 
                "7 Wageningen University, Laboratory of Phytopathology, Wageningen, The Netherlands", br(), 
                "8 Marine Sciences and Applied Biology, University of Alicante, Alicante, Spain", br(), 
                "9 Addis Ababa University, Addis Ababa, Ethiopia", br(), 
                "10 Experimental and Computational Plant Development, Utrecht University, Padualaan 8, 10 3584 CH, The Netherlands", br(), 
                "11 Plant Stress Resilience, Utrecht University, Padualaan 8, 10 3584 CH, The Netherlands", br(), 
                "12 Ecology and Biodiversity, Institute of Environmental Biology, Utrecht University, Padualaan 8, 10 3584 CH, Utrecht, The Netherlands", br(), 
               "13 Howard Hughes Medical Institute, University of California, Davis, Davis, CA 95616", br(), 
                "* Equal contributions"),
            h5("Jiayin Song and Shaina Rahman contributed to early Shiny App code development during the Young Scholar’s Program of 2021."),
            h5("This work was supported, in whole or in part, by the Gates Foundation grant to Netherlands Institute of Ecology (NIOO) “Engineering Soil and Plant Microbiomes for Enhanced Crop Productivity in Africa (OPP1082853) and “PROMISE II: Promoting microbes for integrated Striga eradication” (INV-047255). The conclusions and opinions expressed in this work are those of the author(s) alone and shall not be attributed to the Foundation. "),
            h5(""),
            h4("Please report issues at", tags$a(href="https://github.com/TameraTaylor/EMIC-Ethiopian-Soil-Analysis", "Github"))
    ),
    
   tabItem(tabName = "map",
      fluidRow(
        box(width=8,
            h3("Interactive Map of Ethiopia"),
            plotlyOutput("mainmap", width = "100%", height = "1000px")), #when height=100% is not big enough for legend
        box(title="Choose a trait to display", width=4,
            selectInput("colorDet", "Variable Displayed", choices = choices))
      )#close fluid reminder that width must sum to 12
    ),
    
    tabItem(tabName = "soilcomposition",
      fluidRow(
        box(width=12,
            h3("Sample Composition Explorer"),
            h4("Select a point on the map to display its composition below"),
            plotlyOutput("piechartmap", width = "650px", height = "650px"),
            h4("Currently viewing sample: ", textOutput("e")))
      ),#close fluid
      fluidRow(
        box(title="Physical Composition", width=4,
            plotOutput("pie1", width = "80%", height = "500px"),
            downloadButton("downloadpie1", "Save Physio Chart")),
        box(title="Nutrient Composition", width=4,
                   plotOutput("nutrients1", width = "100%"),
                   plotOutput("nutrients2", width = "100%"),
                   plotOutput("nutrients3", width = "100%"),
                   downloadButton("downloadpie2", "Save Chemical Charts")), #all 3
        box(title="Microbial Composition", width=4,
            plotOutput("pie3", width = "100%"),
            plotOutput("pie4", width = "100%"),
            downloadButton("downloadpie3", "Save Microbe Chart"))
      ) #close fluid
    ), #close soil comp
    
    tabItem(tabName = "correlationanalysis",
      fluidRow(
        column(width=4,
               box(width=NULL, 
                   pickerInput(
                     inputId = "ID", 
                     label = "Choose the samples you want to see:", 
                     choices = data$ID,
                     selected = data$ID, 
                     multiple = T,
                     width = "200px",
                     options = list(
                       `actions-box` = TRUE,
                       `deselect-all-text` = "None",
                       `select-all-text` = "All",
                       `none-selected-text` = "zero"
                     )),
                   virtualSelectInput(
                     inputId = "variable",
                     label = "Choose the properties you want to see:",
                     choices = list(`Local` = list("Altitude","lon","lat","Agro_Climate_type"),
                                    `Striga` = list("Normalized_Striga_Count_Ave","Normalized_Striga_Count_SE",
                                                    "Seedbank_Ave","Seedbank_SE","log.Seedbank_Ave",
                                                    "log.Seedbank_SE","log.Infestation_Ave","log.Infestation_SE"),
                                    `Physical` = list("pH","organic","Carbonated_lime","Clay","Silt","Sand","Clay_humus"),
                                    `Chemical` = list("Mg_avail","Mg_total","Mg",
                                                      "Na_avail","Na_total","Na",
                                                      "K_avail","K_total","K",
                                                      "Ca_avail","Ca_total","Ca",
                                                      "P_avail","P_total","S_avail",
                                                      "S_total","S.delivery","C.S_ratio",
                                                      "C.OS.ratio","C","C_inorganic",
                                                      "N","C.N_ratio","N.delivery",
                                                      "Microbial_activity"),
                                    `Microbial` = list("FUNGUS_ShannonIndex", "BACTERIA_ShannonIndex")),
                     multiple = TRUE
                   ),
                   selectInput(inputId="graph", "Choose a type of Analysis", 
                               choices = c("Correlation Matrix", 
                                           "Pairs Scatter Plots", 
                                           "Principal Component Analysis",
                                           "Heatmap"), multiple=FALSE),
                   uiOutput("CorrTypeSelect"),
                   actionButton("corrGO", "Go!", class="btn-success"),
                   downloadButton("downloadcorrtab", "Download Analysis"))),
        column(width=8,
               box(width=NULL,  
                   plotOutput("plots", height = "700px"),
                   textOutput("stattext")
                   # tableOutput(outputId= "CorText")
                   ))
      )#close fluid
    ),#close correlationanalysis
    
    tabItem(tabName = "modeling",
      fluidRow(
        column(width=4,
               box(width=NULL, 
                   pickerInput(
                     inputId = "ModID", 
                     label = "Choose the samples you want to see:", 
                     choices = data$ID,
                     selected = data$ID, 
                     multiple = T,
                     width = "200px",
                     options = list(
                       `actions-box` = TRUE,
                       `deselect-all-text` = "None",
                       `select-all-text` = "All",
                       `none-selected-text` = "zero"
                     )),
                   virtualSelectInput(
                     inputId = "ModDVariable",
                     label = "Choose the DEPENDENT variable you want to model with:",
                     choices = list(`Local` = list("Altitude","lon","lat","Agro_Climate_type"),
                                    `Striga` = list("Normalized_Striga_Count_Ave","Normalized_Striga_Count_SE",
                                                    "Seedbank_Ave","Seedbank_SE","log.Seedbank_Ave",
                                                    "log.Seedbank_SE","log.Infestation_Ave","log.Infestation_SE"),
                                    `Physical` = list("pH","organic","Carbonated_lime","Clay","Silt","Sand","Clay_humus"),
                                    `Chemical` = list("Mg_avail","Mg_total","Mg",
                                                      "Na_avail","Na_total","Na",
                                                      "K_avail","K_total","K",
                                                      "Ca_avail","Ca_total","Ca",
                                                      "P_avail","P_total","S_avail",
                                                      "S_total","S.delivery","C.S_ratio",
                                                      "C.OS.ratio","C","C_inorganic",
                                                      "N","C.N_ratio","N.delivery",
                                                      "Microbial_activity"),
                                    `Fungal Family` = list("FUNGUS_ShannonIndex", "Aphelidiomycota","Ascomycota",
                                                           "Basidiobolomycota","Basidiomycota","Blastocladiomycota",
                                                           "Calcarisporiellomycota","Chytridiomycota","Entomophthoromycota",
                                                           "Entorrhizomycota","Fungi_phy_Incertae_sedis","Glomeromycota",
                                                           "Kickxellomycota","Monoblepharomycota","Mortierellomycota",
                                                           "Mucoromycota","Neocallimastigomycota","Olpidiomycota",
                                                           "Rozellomycota","Fungus_NA"),
                                    `Bacterial Family`=list("BACTERIA_ShannonIndex","Abditibacteriota",
                                                            "Acidobacteriota","Actinobacteriota","Aenigmarchaeota",
                                                            "Armatimonadota","Bacteroidota","Bdellovibrionota",
                                                            "Chloroflexi","Crenarchaeota","Cyanobacteria",
                                                            "Dadabacteria","Deinococcota","Dependentiae",
                                                            "Desulfobacterota","Elusimicrobiota","Entotheonellaeota",
                                                            "FCPU426","Fibrobacterota","Firmicutes",
                                                            "GAL15","Gemmatimonadota","Hydrogenedentes",
                                                            "Latescibacterota","MBNT15","Methylomirabilota",
                                                            "Myxococcota","Bacteria_NA","NB1.j",
                                                            "Nanoarchaeota","Nitrospirota","Patescibacteria",
                                                            "Planctomycetota","Proteobacteria","RCP2.54",
                                                            "SAR324.clade.Marine.group.B.","Spirochaetota","Sumerlaeota",
                                                            "TX1A.33","Thermoplasmatota","Thermotogota",
                                                            "Verrucomicrobiota","WS2","Zixibacteria")),
                     multiple = FALSE
                   ),
                   virtualSelectInput(
                     inputId = "ModIVariable",
                     label = "Choose the INDEPENDENT variables you want to model with:",
                     choices = list(`Local` = list("Altitude","lon","lat","Agro_Climate_type"),
                                    `Striga` = list("Normalized_Striga_Count_Ave","Normalized_Striga_Count_SE",
                                                    "Seedbank_Ave","Seedbank_SE","log.Seedbank_Ave",
                                                    "log.Seedbank_SE","log.Infestation_Ave","log.Infestation_SE"),
                                    `Physical` = list("pH","organic","Carbonated_lime","Clay","Silt","Sand","Clay_humus"),
                                    `Chemical` = list("Mg_avail","Mg_total","Mg",
                                                      "Na_avail","Na_total","Na",
                                                      "K_avail","K_total","K",
                                                      "Ca_avail","Ca_total","Ca",
                                                      "P_avail","P_total","S_avail",
                                                      "S_total","S.delivery","C.S_ratio",
                                                      "C.OS.ratio","C","C_inorganic",
                                                      "N","C.N_ratio","N.delivery",
                                                      "Microbial_activity"),
                                    `Fungal Family` = list("FUNGUS_ShannonIndex", "Aphelidiomycota","Ascomycota",
                                                           "Basidiobolomycota","Basidiomycota","Blastocladiomycota",
                                                           "Calcarisporiellomycota","Chytridiomycota","Entomophthoromycota",
                                                           "Entorrhizomycota","Fungi_phy_Incertae_sedis","Glomeromycota",
                                                           "Kickxellomycota","Monoblepharomycota","Mortierellomycota",
                                                           "Mucoromycota","Neocallimastigomycota","Olpidiomycota",
                                                           "Rozellomycota","Fungus_NA"),
                                    `Bacterial Family`=list("BACTERIA_ShannonIndex","Abditibacteriota",
                                                            "Acidobacteriota","Actinobacteriota","Aenigmarchaeota",
                                                            "Armatimonadota","Bacteroidota","Bdellovibrionota",
                                                            "Chloroflexi","Crenarchaeota","Cyanobacteria",
                                                            "Dadabacteria","Deinococcota","Dependentiae",
                                                            "Desulfobacterota","Elusimicrobiota","Entotheonellaeota",
                                                            "FCPU426","Fibrobacterota","Firmicutes",
                                                            "GAL15","Gemmatimonadota","Hydrogenedentes",
                                                            "Latescibacterota","MBNT15","Methylomirabilota",
                                                            "Myxococcota","Bacteria_NA","NB1.j",
                                                            "Nanoarchaeota","Nitrospirota","Patescibacteria",
                                                            "Planctomycetota","Proteobacteria","RCP2.54",
                                                            "SAR324.clade.Marine.group.B.","Spirochaetota","Sumerlaeota",
                                                            "TX1A.33","Thermoplasmatota","Thermotogota",
                                                            "Verrucomicrobiota","WS2","Zixibacteria")),
                     multiple = TRUE
                   ),
                   actionButton("modelingGO", "Go!", class="btn-success"))),
        column(width=8,
               box(width=NULL, 
                   #tableOutput("forestdf"),
                   plotOutput("histogram"),
                   plotOutput("tree"),
                   verbatimTextOutput("treestats"),
                   verbatimTextOutput("randomforest"),
                   plotOutput("treeimportance"),
                   textOutput("treetext")))
      )#close fluid
    ), #close modeling
   
   
    tabItem(tabName="Data_Browser", #just display the datatables of ASVs & taxonomy
            h2("Walk through the microbial data preparation"),
        fluidRow(
          column(width=4, selectInput(inputId="Taxa1", label="Choose the kingdom level to view:", choices=c("Bacteria", "Fungi"), multiple = FALSE)),
          column(width=4,selectInput(
            inputId = "Taxa2", 
            label = "Choose at what taxa level you would like to view the data:", 
            choices = c("phylum", "class", "order", "family", "genus", "ASV"),
            multiple = FALSE)), 
          column(width=4, pickerInput(
            inputId = "ID2", 
            label = "Choose the samples you want to see:", 
            choices = list("E01",     "E02",     "E03",     "E04",     "E05",     "E06",     "E07",     "E08",    
                        "E09",     "E10",     "E11",     "E12",     "E13",     "E14",     "E15",     "E16",     "E17",    
                        "E18",     "E19",     "E20",     "E21",     "E22",     "E23",     "E24",     "E25",     "E26",    
                        "E27",     "E28",     "E29",     "E30",     "E31",     "E32",     "E33",     "E34",     "E35",    
                        "E36",     "E37",     "E38",     "E39",     "E40",     "E41",     "E42",     "E43",     "E44",    
                        "E45",    "E46",     "E47",     "E48"),
            selected = list("E01",     "E02",     "E03",     "E04",     "E05",     "E06",     "E07",     "E08",    
                            "E09",     "E10",     "E11",     "E12",     "E13",     "E14",     "E15",     "E16",     "E17",    
                            "E18",     "E19",     "E20",     "E21",     "E22",     "E23",     "E24",     "E25",     "E26",    
                            "E27",     "E28",     "E29",     "E30",     "E31",     "E32",     "E33",     "E34",     "E35",    
                            "E36",     "E37",     "E38",     "E39",     "E40",     "E41",     "E42",     "E43",     "E44",    
                            "E45",    "E46",     "E47",     "E48"),
            multiple = T,
            # width = "200px",
            options = list(
              `actions-box` = TRUE,
              `deselect-all-text` = "None",
              `select-all-text` = "All",
              `none-selected-text` = "zero") ) 
            ) #close column
          ), #close fluidrow
        fluidRow(actionButton("MicrobeButton1", "Confirm", class="btn-success")),
        fluidRow(
          plotOutput(outputId= "TaxaAbundance"),
          downloadButton("MicrobeBarDownload", "Download Plot") ),
        fluidRow(
          h4("Preview of the raw data summed up at your chosen taxa level."),
          div(style = 'overflow-x: scroll', DT::dataTableOutput("RawChosenTaxa")),
          downloadButton("RawChosenTaxaDownload", "Download Full Data") ),
        fluidRow(
          h4("Preview of the data after imputation (zeros were replaced with  a random number between 0 and the detection limit) and centered log-ratio transformation."),
          uiOutput("CLRChosenTaxa"),
          downloadButton("CLRChosenTaxaDownload", "Download Full Data") )         
     ),
        
        
     tabItem(tabName="Microbe_Analysis",
        fluidRow(
          column(width=4,
                 box(width=NULL,
                   textOutput("MicrobeChoices")
                 )
          ),
          column(width=4,
                 box(width=NULL,
                     virtualSelectInput(
                       inputId = "variable2",
                       label = "Choose the properties to correlate with your prepared microbe data:",
                       choices = list(`Local` = list("Altitude","lon","lat","Agro_Climate_type"),
                                      `Striga` = list("Normalized_Striga_Count_Ave","Normalized_Striga_Count_SE",
                                                      "Seedbank_Ave","Seedbank_SE","log.Seedbank_Ave",
                                                      "log.Seedbank_SE","log.Infestation_Ave","log.Infestation_SE"),
                                      `Physical` = list("pH","organic","Carbonated_lime","Clay","Silt","Sand","Clay_humus"),
                                      `Chemical` = list("Mg_avail","Mg_total","Mg",
                                                        "Na_avail","Na_total","Na",
                                                        "K_avail","K_total","K",
                                                        "Ca_avail","Ca_total","Ca",
                                                        "P_avail","P_total","S_avail",
                                                        "S_total","S.delivery","C.S_ratio",
                                                        "C.OS.ratio","C","C_inorganic",
                                                        "N","C.N_ratio","N.delivery",
                                                        "Microbial_activity"),
                                      `Microbial` = list("FUNGUS_ShannonIndex", "BACTERIA_ShannonIndex")),
                       multiple = TRUE)
                      ) ),
          column(width=4,
                 box(width=NULL,
                     selectInput(inputId="CorrType2", "Choose correlation type:",
                                 choices=c("pearson", "spearman"), multiple=FALSE),
                     uiOutput("CorrTypeSelect2"),
                     actionButton("MicrobeButton2", "Go!", class="btn-success")
                     ) )
        ), #fluidrow
        fluidRow(
          # div(style='overflow-x: scroll', DT::dataTableOutput("microbeAnalysis")),
          h3("Correlation Coefficients of your analysis."),
          div(style='overflow-x: scroll', DT::dataTableOutput("microbeCorr")),
          h4("Click to download the excel workbook of these analysis results."),
          downloadButton("downloadcorrtab2", "Download Analysis")
        ),#fluidrow
        fluidRow(
            title = "Correlation Plot", width = 12, solidHeader = TRUE, status = "primary",
            plotOutput("MicrobeCorrelogram", height = "600px"),
            downloadButton("downloadMicroCorrPlot", "Download Correlation Plot")
        )
        
      ) #tab
  )#close tabItems
)#close body

# Put them together into a dashboardPage
ui<- dashboardPage(skin = "purple",
  dashboardHeader(title = "PROMISE Analysis App"),
  sidebar,
  body)

