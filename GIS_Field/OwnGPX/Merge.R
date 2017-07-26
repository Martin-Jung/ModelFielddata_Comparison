#### Preperation and standard Package Loading ####
rm(list=ls())
library(rgdal)
library(gdata)
l <- list.files(".","*.gpx")
asP <- function(x) readOGR(x,layer="waypoints",verbose=F)
numbers_from_string <- function(x) as.numeric(gsub("\\D", "", x))
out <- lapply(l,asP)
# Loop through all the gpx files, subset the right columns 
res <- SpatialPointsDataFrame(out[[1]],as.data.frame(out[[1]]))
res$PlotID <- numbers_from_string(res$name)
res <- subset(res,select=c("ele","time","name","PlotID"))
for(i in 2:length(out)){
  w <- out[[i]]
  w$PlotID <- numbers_from_string(w$name)
  s <- subset(w,select=c("ele","time","name","PlotID"))
  res <- rbind(res,s)
  rm(w,s)
}
res <- res[which(res$name!="KILI KIDIA"),];res <- res[which(res$name!="STATION"),]
# Kickout Sample T006 as it was merged with Plot01 due to Proximity
res <- res[which(res$name!="006"),]
res <- res[which(res$name!="001"),] # And kick out test-coordinate from nairobi

# Make correct PLOTID
res$name <- as.character(res$name)
res$name[which(startsWith(res$name,pattern="T"))] <- "001"
res$PlotID <- paste0(ifelse(res$PlotID > 053,"K","T"),res$name)

# Deal with the double name by adding a b to it
res$PlotID[which(res$name=="082" & res$ele < 1000)] <- paste0(res$PlotID[which(res$name=="082" & res$ele < 1000)],"b")
res$name[which(res$name=="082" & res$ele < 1000)] <- paste0(res$name[which(res$name=="082" & res$ele < 1000)],"b")

# Add my name
res$Authority <- "MartinJ"

# Load in Dickens plots and add them as well
dick <- readOGR("../dickens_plots/","Dickens_Merged")
names(dick) <- c("ele","PlotID")
dick$Authority <- "DickensO"
dick$time <- NA; dick$name <- NA
out <- rbind(res,dick)

writeOGR(out,dsn=".",layer="PlotLocalities",driver="ESRI Shapefile",overwrite_layer=T,verbose=F)
