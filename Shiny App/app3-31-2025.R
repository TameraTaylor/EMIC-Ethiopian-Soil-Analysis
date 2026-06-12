#Required packages
pacman::p_load(shiny, shinydashboard, shinythemes, shinyWidgets, #shiny app
           reshape2, tidyr, dplyr, stringr, data.table, #data wrangling 
           ggplot2, ggspatial, ggfortify, ggiraph, plotly, #plotting
           stars, # spatiotemporal data handling
           raster, # raster data handling
           terra, # raster data handling
           sf, # vector data handling
           geojsonio, "mapproj", rpart.plot, rpart, sp,
           lubridate, # dates handling
           patchwork, # arranging figures
           tigris, # county border
           colorspace, # color scale
           viridis, # arranging figures
           corrplot, "psych", cocor,
           broom, janitor, RColorBrewer, gridExtra,
           pheatmap, vroom, 
           randomForest,
           EBImage
)


###REQUIRED TO UPDATE APP ON SHINYAPPS,IO
#library(rsconnect)
#options(repos = BiocManager::repositories())
#rsconnect::deployApp("C:/Users/small/Documents/Brady lab/Shiny App/WebHosting")
###################

########################

#Global functions/actions

#####################################

data <- read.csv("SoilData.csv")
BacteriaData<-read.csv("BACTERIA ASVs.csv", row.names = "ASV_ID")
FungiData<-read.csv("FUNGUS ASVs.csv", row.names = "ASV_ID")
SoilData_units <- read.csv("SoilData_units.csv")




#________________________________MAP_____________
eth <- geojson_read("https://raw.githubusercontent.com/georgique/world-geojson/master/countries/ethiopia.json", what = "sp")
AgroEcors <- stack("KG3_CRUTS32_Hist_8110raster.tif")

## crop and mask for just ethiopia shape
r2 <- crop(AgroEcors, extent(eth))
r3 <- mask(r2, eth)

ethiopia_rds <- file.path("gadm36_ETH_3_sf.rds")
ethiopia_sf <- readRDS(ethiopia_rds)
ethiopia_regions_sf <-
  ethiopia_sf %>%
  mutate(NAME_1 = gsub("Southern Nations, Nationalities and Peoples", "SNNP",
                       NAME_1)) %>%
  group_by(NAME_1) %>%
  summarise() %>%
  ungroup() %>%
  st_as_sf()

#Getting the colors
r4 <-as.data.frame(r3, xy = TRUE) %>%
  na.omit() #--- remove cells with NA for any of the layers ---#
r4$KG3_CRUTS32_Hist_8110raster<-as.factor(r4$KG3_CRUTS32_Hist_8110raster)

mypalette<-colorRampPalette(rev(brewer.pal(9,"YlGn")))(14)

#Base Ethiopia Map
MapBase<-  ggplot(data=r4) +
  geom_raster(aes(x = x, y = y, fill = KG3_CRUTS32_Hist_8110raster)) +
  geom_sf(data=ethiopia_regions_sf, fill=NA, linewidth=.8)+
  scale_fill_manual(values=rev(mypalette))+ 
  labs(fill="Agroclimate Zone")+
  theme_void() 
#____________________




choices <- c()
for(i in 1:length(data)){
  if(is.numeric(data[,i]))
    choices <- c(choices, colnames(data[i]))
}








#######################################RUN APP

source('ui 3-31-2025.R', local = TRUE)
source('server3-31-2025.R', local=TRUE)
shinyApp(ui, server)
