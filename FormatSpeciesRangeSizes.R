# Short script to format species range sizes
library(rgdal)
library(sp)
library(raster)
library(parallel)
library(rgeos)
library(dplyr)
library(GGally)

allsp_africa <- readOGR("GIS","Spp_1DegCells_clip")
atiles <- readOGR("GIS","GBIF_Africa1degCells")
d <- as.data.frame(allsp_africa)

# Calculate range based on global Shape_Area
res <- d %>% dplyr::select(SCINAME,Shape_Area) 

# Calculate area based on mollweide equal area projection
sub <- spTransform(allsp_africa,CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")) 
res$Shape_Area_Africa <- gArea(sub,byid=T)


# Calculate Area based on occupied Gridcells
res$DegOccupiedArea <- NA
for(species in unique(allsp_africa$SCINAME)){
  print(species)
  sub <- subset(allsp_africa,SCINAME==species)
  a <-  which(gCoveredBy(atiles,sub,byid=T))
  res$DegOccupiedArea[which(res$SCINAME==species)] <- length(a)
}

# Load the GBIF Data
a <- readRDS("GBIF_OCC-Fielddata.rds")
b <- readRDS("GBIF_OCC-PREDICTS.rds")
d <- rbind(a,b) %>% select(SpeciesID,NrCells)
d[-which(duplicated(d)),]

ggpairs(res,3:ncol(res),
        upper = list(continuous = "cor", combo = "dot"),
        lower = list(continuous = "points", combo = "dot"))


saveRDS(res,"BirdLifeRangeSizes.rds")

