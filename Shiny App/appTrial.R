library(shiny)
library(shinydashboard)
library(shinythemes)
library(shinyWidgets)
library(DT)
library(reshape2)
library(tidyr)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(ggspatial)
library(ggfortify)
library(ggiraph)
library(ggrepel)
library(plotly)
library(stars)
library(raster)
library(terra)
library(sf)
library(geojsonio)
library(mapproj)
library(rpart.plot)
library(rpart)
library(sp)
library(lubridate)
library(patchwork)
library(tigris) 
library(colorspace)
library(viridis)
library(corrplot)
library(psych)
library(cocor)
library(broom)
library(janitor)
library(RColorBrewer)
library(gridExtra)
library(pheatmap)
library(vroom)
library(randomForest)
library(EBImage)
library("Tjazi")
library(openxlsx)




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

#microbe graphs
mycolors<-c(
  "#faf600",
  "#f09341",
  "#f04141",
  "#f041aa",
  "#c741f0",
  "#737ffa",
  "#41c1f0",
  "#00fad5",
  "#2cc946",
  "#abfa00",
  "#A3A3A3")

#######################################RUN APP

source('scripts/ui.R', local = TRUE)
source('scripts/server.R', local=TRUE)
shinyApp(ui, server)
