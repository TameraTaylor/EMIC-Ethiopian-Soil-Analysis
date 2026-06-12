#separate ui for shiny app
print(search())
print(sessionInfo()$otherPkgs)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Home", tabName = "home"),
    menuItem("Ethiopian Soil Summary",
             menuSubItem("Map", tabName = "map"),
             menuSubItem("Soil Composition", tabName = "soilcomposition"),
             menuSubItem("Correlation Analysis", tabName = "correlationanalysis"),
             menuSubItem("Modeling", tabName = "modeling")),
    menuItem("Ethiopian Soil Microbiome",
             menuSubItem("Data Browser"), #just display the datatables of ASVs & taxonomy
             menuSubItem("Correlation Analysis")
             )
    )
  )




body <- dashboardBody(
  tabItems(
    tabItem(tabName = "home", 
            h2("Exploration of Ethiopian Soil Data from PROMISE"),
            h3("PROMISE (Promoting Root Microbes against Striga in Ethiopia) has collected these data as a collaborative effort between institutions.",
               "This application has been designed to offer means of analysis of these data presented in ________ to a wider audience.",
               "Future plans are to adapt the application to analyze user inputted data and provide treatment reccomdations for striga infested soils."),
            imageOutput("photo1"), #also add logos of collaborators
            h5("Representative photo of sorghum field in Ethiopia sampled by PROMISE. Photo courtesy: Dr. Taye Tessema"),
            h2("Description"),
            h4("A total of 48 composite soil samples covering a trajectory of more than 1500 km were collected from naturally Striga infested sorghum growing agro ecological zones in Ethiopia.",
               "Soil samples were collected from naturally Striga infested sorghum fields in Amhara (Kemise, North Shewa, South and North Wollo Zones) and Tigray (West, Central and South zones) regions of Ethiopia.",
               "Sorghum fields with four categories (zero, low, medium and high) of Striga field infestation were randomly selected, and samples were taken from the top layer (0-20 cm) around the root zone of ",
               "a sorghum plant then combined together to form one composite sample per field.",
               "Physical, chemical, and microbiological compositions were analyzed for these samples in a joint effort by the PROMISE Project."),
            h4("[INSERT READ ME TEXT HERE]")
    ),
    
    tabItem(tabName = "map",
      fluidRow(
        box(title="Map of Ethiopia", width=8,
            plotlyOutput("mainmap", width = "100%", height = "1000px")), #when height=100% is not big enough for legend
        box(title="Attribute Display", width=4,
            selectInput("colorDet", "Variable Displayed", choices = choices))
      )#close fluid reminder that width must sum to 12
    ),
    
    tabItem(tabName = "soilcomposition",
      fluidRow(
        box(title="Map of Ethiopia", width=12,
            h4("Select a point on the map to display its composition"),
            plotlyOutput("piechartmap", width = "650px", height = "650px"),
            h4(textOutput("e")))
      ),#close fluid
      fluidRow(
        box(title="Physical Composition", width=4,
            plotOutput("pie1", width = "80%", height = "500px"),
            downloadButton("downloadpie1", "Save Physio Chart")),
        box(title="Nutrient Composition", width=4,
            column(width = 4, plotOutput("nutrients1", width = "100%"),
                   plotOutput("nutrients2", width = "100%"),
                   plotOutput("nutrients3", width = "100%"),
                   downloadButton("downloadpie2", "Save Chemical Charts"))), #all 3
        box(title="Microbial Composition", width=4,
            plotOutput("pie3", width = "100%", height = "500px"),
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
                   textOutput("stattext"),
                   tableOutput(outputId= "CorText")))
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
    ) #close modeling
        #REMOVED THE SORGHUM DESCRIPTION THINGS
  )#close tabItems
)#close body

# Put them together into a dashboardPage
ui<- dashboardPage(skin = "purple",
  dashboardHeader(title = "PROMISE Analysis App"),
  sidebar,
  body)

