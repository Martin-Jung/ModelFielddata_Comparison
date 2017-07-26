#### Cleaning
par.ori <- par(no.readonly=T)
library(gdata)
library(raster)
library(rgdal)
library(lubridate)
library(plyr)
library(stringr)
library(MODISTools)
library(maptools)
library(MODIS)  
ori.dir <- "/home/martin/Documents/Studium/Kopenhagen/Masters/Analysis_PREDICTS"
setwd(ori.dir)
pd.data.dir <- "Data/"
fd.data.dir <- "../FieldData/"
gis.folder <- "../GIS_Field/"
gis.folder2 <- "../GIS/"
gis.extern <- "/media/GISbox/Raster/"
#fff
field_data_prep <- function(sites){
  stopifnot(file.exists(gis.extern))
  # Make Results Dataframe
  stop("Was initialy executed with wrong coordinates. Lat-long switched")
  d <- data.frame(ID = sites$ID,pop=NA,elev=NA,meanEVI=NA,yieldEVI=NA,meanNDVI=NA,yieldNDVI=NA,FC2000=NA)
  
  # Prepare Spatial Data
  #t_buff <- readOGR(gis.folder ,"Taita_transectbuffer_1km_Arc1960_37s",verbose=F)
  #k_buff <- readOGR(gis.folder ,"Kilimanjaro_transect_buffer1km_arc1960",verbose=F)
  #kili_lulc <- raster(paste0(gis.folder,"CHIESA_Landuse/lulc_clutter_kili_transect_2012_final.img"))
  
  # ---------------- #
  #### Load DEMs
  # DEM not needed as Data directly quiered from GPS
  # But for sake of conistency we sample from the same source as the large-scale data (PREDICTS)
  # Get Elevation data for all sites SOURCE = SRTM
  # Instead of Downloading we query and untar the respective grids from 
  # the external harddrive
  tiles <- readOGR(paste0(gis.extern,"/DEM/SRTM/"),"srtm-5dg-grid")
  srtm <- paste0(gis.extern,"DEM/SRTM/srtm.csi.cgiar.org/SRTM_v41/SRTM_Data_GeoTIFF/")
  sp <- SpatialPointsDataFrame(cbind(sites$Long,sites$Lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),data=sites)
  ov <- over(sp,y=tiles);ov$file <- paste0(ov$NAME,".zip") # Query SRTM Tile
  # Loop through all respective files on the harddrive,
  # extract the values and append them to pred_gisdata  
  for(id in 1:nrow(ov)){
    print(paste0("Processing NR:",id," and tile ",ov$file[id]))
    f = paste0(srtm,ov$file[id])
    ras = paste0("/tmp/SRTM/",ov$NAME[id],".tif")
    # Clear up Temp if different and unzip file
    if(!file.exists(ras)){      
      unlink("/tmp/SRTM/*")     
      unzip(zipfile=f,exdir="/tmp/SRTM",overwrite=T)
      dem = raster(ras)      
    }
    d$elev[id] <- raster::extract(dem,coordinates(sp[id,]))
    rm(ras)    
    print(paste0(nrow(ov)-id," remaining..."))    
  }  
  rm(tiles,srtm,ov,sp,id,dem,f);unlink("/tmp/SRTM/*")
  
  # ------------------#
  #### Population
  popdir <- paste0(gis.extern,"Anthropological/Population/Detailed/")  
  sp <- SpatialPointsDataFrame(cbind(sites$Long,sites$Lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),data=sites)
  data(wrld_simpl)
  #first identify failures
  country <- ifelse(sites$Transect=="Taita","Kenya","United Republic of Tanzania")
  matchLookup <- match(as.character(country),as.character(wrld_simpl$NAME))
  # Unzip all files - messy - big - takes long
  dir.create("/media/GISbox/Temp") # Create on temp because temporary is to small
  unlink("/media/GISbox/Temp/*")
  for(id in 1:length(matchLookup)){
    code = wrld_simpl$ISO3[matchLookup[id]]
    f = paste0(popdir,code,"-POP.7z")
    if(!file.exists(f)) cat(code," does not exist as file on the drive")
    if( length( list.files("/media/GISbox/Temp/",pattern=as.character(code)) )==0){
      cat("Extracting ",code,"\n")
      system(paste0("7z e -y -o/media/GISbox/Temp/ ",f," > /dev/null"))
      cat("Finished extracting of file ",f,"\n")
    }
  }
  
  # Kick out all that are not pure tifs
  unlink(list.files("/media/GISbox/Temp/",pattern="*.tif.",full.names=T))
  unlink(list.files("/media/GISbox/Temp/",pattern=".tfw",full.names=T))
  # Kickout 2015 projections and other projections
  unlink(grep(x=list.files("/media/GISbox/Temp/",pattern="*.tif",full.names=T),pattern="10",value=T,invert=T))
  # Manually delete dublicated files
  unlink(list("/media/GISbox/Temp/UGA_pph_v2b_2010_UNadj.tif","/media/GISbox/Temp/NAM_pph_v2b_2010_UNadj.tif"))
  ff=0
  for(id in 1:length(matchLookup)){
    code = wrld_simpl$ISO3[matchLookup[id]]
    f = grep(x=list.files("/media/GISbox/Temp/",pattern="adj",full.names=T),pattern=code,value=T)
    if(length(f)>1){
      stop(paste("There are more than one adj file for", code))      
    } else {
      if(f != ff)  pop = raster(f)
      ff = f
      ex <- try( raster::extract(pop,coordinates(sp[id,])) )
      if (class(ex)[1]!= "try-error") {
        d$pop[id] <- ex        
        print(paste0(length(matchLookup)-id," remaining..."))
      } else {
        print(paste0(id," could not be queried. Raster missing..."))
      }  
    } 
  }  
  rm(f,id,code,pop,wrld_simpl)
  unlink("/media/GISbox/Temp/*");file.remove("/media/GISbox/Temp")
  # Taita didn't work for some reason. Using source from computer
  #t.pop <- raster("../GIS/Afripop/ke10_crop.tif")
  #df$pop[which(is.na(df$pop))] <- raster::extract(t.pop,sp[which(is.na(df$pop)),])
  #rm(t.pop)
  
  
  # ----------- #
  #### MODIS ####
  # MODISTools with pixelwise past 2 years as well  
  cl_sub <- sites
  # MODISTools requires you to make a query data.frame
  coord <- cbind(cl_sub$Lat,cl_sub$Long)
  dates <- cbind(year(cl_sub$Sample_start_earliest),year(cl_sub$Sample_end_latest))
  dates[,1] <- dates[,1] - 2
  dates[,2] <- ifelse(dates[,1]<2000,2002,dates[,2])
  dates[,1] <- ifelse(dates[,1]<2000,2000,dates[,1])

  product <- "MOD13Q1"
  bands <- c("250m_16_days_NDVI","250m_16_days_EVI","250m_16_days_pixel_reliability") # What to query. You can get the names via GetBands
  savedir <- "PredictsMODIS2years_FW/" # You can save the downloaded File in a specific folder
  if(!file.exists(savedir)) dir.create(savedir)
  pixel <- c(0,0) # Get the central pixel only (0,0) or a quadratic tile around it
  #which(!duplicated(sites.pred$Source_ID))
  i <- 1
  period <- data.frame(lat=coord[i:nrow(coord),1],long=coord[i:nrow(coord),2],start.date=dates[i:nrow(coord),1],end.date=dates[i:nrow(coord),2],id=cl_sub$ID)
  period <- UpdateSubsets(period,StartDate = T,Dir = savedir)
  # To download the pixels  
  k = 1
  MODISSubsets(LoadDat = period[k:nrow(period),],Products = product,Bands = bands,Size = pixel,SaveDir = savedir,StartDate = T)
  
  # Load summary file
  #source("MODIS_YieldCalculation.R")
  #summaries <- MODIS_yield(savedir,period,c(0))
  #summaries$missing
  MODISSummaries(LoadDat = period, Product = "MOD13Q1",Dir = savedir,
                 Bands = c("250m_16_days_EVI","250m_16_days_NDVI"),
                 ValidRange = c(-2000,10000), NoDataFill = -3000, ScaleFactor = 0.0001,
                 StartDate = TRUE, QualityScreen = TRUE, QualityThreshold = 0,
                 QualityBand = "250m_16_days_pixel_reliability",
                 Interpolate=T,Mean=T,Yield=T)
  
  mr <- read.csv("PredictsMODIS2years_FW/MODIS_Summary_MOD13Q1_2014-10-01_h17-m2-s36.csv",header=T)
  md <- read.csv("PredictsMODIS2years_FW/MODIS_Data_MOD13Q1_2014-10-01_h17-m2-s36.csv",header=T)
  stopifnot(md$id==period$id)
  stopifnot(length(which(!(period$id%in%md$id)))==0)
  
  # CAREFUL! Longitude / Latitude switched for some samples
  d <- readRDS("Data/SiteData.rds")
  m_ev <- subset(mr,data.band=="250m_16_days_EVI")
  m_ev$start.date <- year(m_ev$start.date)
  m_ev$end.date <- year(m_ev$end.date)
  m_ev$index <- paste(round(m_ev$long,5),round(m_ev$lat,5),m_ev$start.date,m_ev$end.date,sep="_")
  period$index <- paste(round(period$long,5),round(period$lat,5),period$start.date,period$end.date,sep="_")
  m <-join(period,m_ev,by="index")
  
  #d$meanEVI[which(!is.na(m$mean.band))] <- m$mean.band[which(!is.na(m$mean.band))]
  #d$yieldEVI[which(!is.na(m$mean.band))] <- m$band.yield[which(!is.na(m$mean.band))]
  d$meanEVI <- m$mean.band
  d$yieldEVI <- m$band.yield
  
  m_ev <- subset(mr,data.band=="250m_16_days_NDVI")
  m_ev$start.date <- year(m_ev$start.date)
  m_ev$end.date <- year(m_ev$end.date)
  m_ev$index <- paste(round(m_ev$long,5),round(m_ev$lat,5),m_ev$start.date,m_ev$end.date,sep="_")
  period$index <- paste(round(period$long,5),round(period$lat,5),period$start.date,period$end.date,sep="_")
  m <-join(period,m_ev,by="index",match = "first")
  
  #d$meanNDVI[which(!is.na(m$mean.band))] <- m$mean.band[which(!is.na(m$mean.band))]
  #d$yieldNDVI[which(!is.na(m$mean.band))] <- m$band.yield[which(!is.na(m$mean.band))]
  d$meanNDVI <- m$mean.band
  d$yieldNDVI <- m$band.yield
  
  
  
  # ---------------- #
  #### FC Hansen ####
  fc2 <- raster("/media/GISbox/Raster/Landcover/HANSEN_2013/Treecover2000/fc2000.vrt")
  d$FC2000 <- raster::extract(fc2,sp)
  d$FC2000[which(is.na(d$FC2000))] <- 0 #replace nodata with zero
  rm(fc2)
  
  # ---------------- #
  #### BIOCLIM baseline data #### 
  # USE Platts data from PREDICTS_wide
  bcd <- as.data.frame(raster::extract(bioclim,sp))
  d <- cbind(d,bcd)
  unlink("/media/GISbox/Temp/*")
  file.remove("/media/GISbox/Temp")
    
  # Save everything
  saveRDS(d,"Data/SiteData.rds")
}

field_data_load <- function(sites){
  # Load prepared data 
  d <- readRDS("Data/SiteData.rds")
  # Local factors - Weather score
  # Contains a basic index for mean weather quality
  match_w <- data.frame(w=unique(sites$WeatherSecond),ID=c(1,2,4,3))
  w <- data.frame(a=match_w$ID[match(sites$WeatherFirst,match_w$w)],b=match_w$ID[match(sites$WeatherSecond,match_w$w)])
  sites$WeatherScore <- rowMeans(w)
  
  # Convert factors to indices using a match up table
  match_lu <- data.frame(LU=unique(sites$PREDICTS.LU),ID=seq(1:5))
  match_lui <- data.frame(LUI=unique(sites$PREDICTS.LUI),ID=seq(1:3))
  
  sites$LU <- match_lu$ID[match(sites$PREDICTS.LU,match_lu$LU)]
  sites$LUI <- match_lui$ID[match(sites$PREDICTS.LUI,match_lui$LUI)]
  rm(match_w,match_lu,w,match_lui)  
  
#   # Finally append predictors to sitedata
#   if("Plantae" %in% sites$Grouping){
#     ind <- which(sites$Grouping=="Aves"&sites$Authority=="DickensO")
#     stopifnot(sites$SiteName[ind]==sites$SiteName[which(sites$Grouping=="Plantae")])
#   # Just append the rows again to the bottom
#     d <- rbind(d,d[ind,])
#     cb <- cbind(sites,d)  
#   } else {
    cb <- cbind(sites,d)
#  }
  
  # Correct yield for field sites
  print("Appended Data. Now correcting NDVI for seasonality.")  
  cb$Latitude <- cb$Lat
  cb$Longitude <- cb$Long
  sp <- subset(cb,select=c("Latitude","Longitude","yieldNDVI","PREDICTS.LU","PREDICTS.LUI"))
  sp <- SpatialPointsDataFrame(cbind(sp$Longitude,sp$Latitude),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),data=sp)
  cb$yield.ndvi.corr <- correctNDVIyield(sp,F)
  rm(sp)
  return(cb)
}


predicts_data_prep <- function(sites.pred){
  stopifnot(file.exists(gis.extern))
  # Make Results Dataframe
  pred_gisdata <- data.frame(ID = sites.pred$Source_ID,IDs=sites.pred$SSS,pop=NA,elev=NA,meanEVI=NA,yieldEVI=NA,meanNDVI=NA,yieldNDVI=NA,FC2000=NA)
  
  # Get Elevation data for all sites SOURCE = SRTM
  # Instead of Downloading we query and untar the respective grids from 
  # the external harddrive
  setwd("/media/GISbox/Raster/DEM/SRTM/")
  tiles <- readOGR(".","srtm-5dg-grid")
  srtm <- "srtm.csi.cgiar.org/SRTM_v41/SRTM_Data_GeoTIFF/"
  sp <- SpatialPointsDataFrame(cbind(sites.pred$Longitude,sites.pred$Latitude),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),data=pred_gisdata)
  ov <- over(sp,y=tiles);ov$file <- paste0(ov$NAME,".zip")
  # Loop through all respective files on the harddrive,
  # extract the values and append them to pred_gisdata
  for(id in 1:nrow(ov)){
    print(paste0("Processing NR:",id," and tile ",ov$file[id]))
    f = paste0(srtm,ov$file[id])
    ras = paste0("/tmp/SRTM/",ov$NAME[id],".tif")
    # Clear up Temp if different and unzip file
    if(!file.exists(ras)){      
      unlink("/tmp/SRTM/*")     
      try( unzip(zipfile=f,exdir="/tmp/SRTM",overwrite=T) )
      try(expr = dem <- raster(ras))
    }
    try( pred_gisdata$elev[id] <- raster::extract(dem,coordinates(sp[id,])) )
    rm(ras)    
    print(paste0(nrow(ov)-id," remaining..."))    
  }  
  rm(tiles,srtm,ov,sp,id,dem,f)
  setwd(ori.dir)
  
  # Get Population density data
  # Worldpop country-wise 100m resolution
  popdir <- ("/media/GISbox/Raster/Anthropological/Population/Detailed/")  
  sp <- SpatialPointsDataFrame(cbind(sites.pred$Longitude,sites.pred$Latitude),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),data=pred_gisdata)
  data(wrld_simpl)
  #first identify failures
  matchLookup <- match(as.character(sites.pred$Country),as.character(wrld_simpl$NAME))
  failedCodes <- sites.pred$Country[is.na(matchLookup)]
  cat(length(failedCodes),"countries failed to match for Population data")
  # Unzip all files - messy - big - takes long
  dir.create("/media/GISbox/Temp") # Create on temp because temporary is to small
  unlink("/media/GISbox/Temp/*")
  for(id in 1:length(matchLookup)){
    code = wrld_simpl$ISO3[matchLookup[id]]
    f = paste0(popdir,code,"-POP.7z")
    if(!file.exists(f)) cat(code," does not exist as file on the drive")
    if( length( list.files("/media/GISbox/Temp/",pattern=as.character(code)) )==0){
      cat("Extracting ",code,"\n")
      system(paste0("7z e -y -o/media/GISbox/Temp/ ",f," > /dev/null"))
      cat("Finished extracting of file ",f,"\n")
  }
  }
  
  # Kick out all that are not pure tifs
  unlink(list.files("/media/GISbox/Temp/",pattern="*.tif.",full.names=T))
  unlink(list.files("/media/GISbox/Temp/",pattern=".tfw",full.names=T))
  # Kickout 2015 projections and other projections
  unlink(grep(x=list.files("/media/GISbox/Temp/",pattern="*.tif",full.names=T),pattern="10",value=T,invert=T))
  # Manually delete dublicated files
  unlink(list("/media/GISbox/Temp//UGA_pph_v2b_2010_UNadj.tif","/media/GISbox/Temp//NAM_pph_v2b_2010_UNadj.tif"))
  
  # Check if everything is alright
  stop("Control if there is only one adj file of population present")
  list.files("/media/GISbox/Temp/",pattern="*.tif",full.names=T)
  
  ff = F
  for(id in 1:length(matchLookup)){
    code = wrld_simpl$ISO3[matchLookup[id]]
    f = grep(x=list.files("/media/GISbox/Temp/",pattern="adj",full.names=T),pattern=code,value=T)
    if(length(f)==1){      
      if(f != ff)  pop = raster(f)
      ff = f
      ex <- try( raster::extract(pop,coordinates(sp[id,])) )
      if (class(ex)[1]!= "try-error") {
        pred_gisdata$pop[id] <- ex        
        print(paste0(length(matchLookup)-id," remaining..."))
      } else {
        print(paste0(id," could not be queried. Raster missing..."))
      }  
    } else {
      print(paste("There are more or none adj file for", code))      
      } 
    } 

  rm(f,id,code,pop,wrld_simpl)
  unlink("/media/GISbox/Temp/*");file.remove("/media/GISbox/Temp")
  
  # MODIS EVI + NDVI
  # Load point-wise annual data.
  #modis = "/media/GISbox/Raster/VegetationIndex/MODIS_Timeseries/africagrids.net/MOD13Q1_250m"
  #sp <- SpatialPointsDataFrame(cbind(pred_gisdata$long,pred_gisdata$lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),data=pred_gisdata)
  #sp <- SpatialPointsDataFrame(cbind(sites$Long,sites$Lat),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),data=sites)

  library("R.utils")
  #dir.create("/media/GISbox/Temp") # Create on temp because temporary is to small
  #unlink("/media/GISbox/Temp/*")
  #pred_gisdata$modEVI <- NA
  #pred_gisdata$modNDVI <- NA
#   for(id in 1:nrow(sp)){
#     start = sp$start.date[id]
#     end = sp$end.date[id]
#     
#       fe = paste0(modis,"/EVI/Annual/","M13EVIA",start,".tif.gz")
#       fn = paste0(modis,"/NDVI/Annual/","M13NDVIA",start,".tif.gz")
#       gunzip(filename=fe,destname=paste0("/media/GISbox/Temp/EVI_",start,".tif"),skip=T,remove=F)
#       cat("Extracted EVI of year",start,"\n")
#       gunzip(filename=fn,destname=paste0("/media/GISbox/Temp/NDVI_",start,".tif"),skip=T,remove=F)
#       cat("Extracted NDVI of year",start,"\n")
#       
#       rasE <- raster(paste0("/media/GISbox/Temp/EVI_",start,".tif"))
#       rasN <- raster(paste0("/media/GISbox/Temp/NDVI_",start,".tif"))
#       cat("Loaded in Raster layers\n")
#       pred_gisdata$modEVI[id] <- raster::extract(rasE,coordinates(sp[id,]))
#       pred_gisdata$modNDVI[id] <- raster::extract(rasN,coordinates(sp[id,]))
#     
#     rm(fe,fn,rasE,rasN)
#     cat("Cleared Cache. New round.")
#   }
#   rm(modis)
#   unlink("/media/GISbox/Temp/*");file.remove("/media/GISbox/Temp")
#   
  ## <--- MODIS Tools ---> ##
  #### Using the MODISTools Package ####
  # MODISTools with pixelwise past 2 years as well  
  cl_sub <- sites.pred
  # MODISTools requires you to make a query data.frame
  coord <- cbind(cl_sub$Latitude,cl_sub$Longitude)
  dates <- cbind(year(cl_sub$Sample_start_earliest),year(cl_sub$Sample_end_latest))
  dates[,1] <- dates[,1] - 2
  dates[,2] <- ifelse(dates[,1]<2000,2002,dates[,2])
  dates[,1] <- ifelse(dates[,1]<2000,2000,dates[,1])
  product <- "MOD13Q1"
  bands <- c("250m_16_days_NDVI","250m_16_days_EVI","250m_16_days_pixel_reliability") # What to query. You can get the names via GetBands
  savedir <- "PredictsMODIS2years_PREDICTS/" # You can save the downloaded File in a specific folder
  if(!file.exists(savedir)) dir.create(savedir)
  stopifnot(dates[,2]-dates[,-1]==0)
  pixel <- c(0,0) # Get the central pixel only (0,0) or a quadratic tile around it
  i <- 1
  period <- data.frame(lat=coord[i:nrow(coord),1],long=coord[i:nrow(coord),2],start.date=dates[i:nrow(coord),1],end.date=dates[i:nrow(coord),2],id=sites.pred$Source_ID[i:nrow(sites.pred)])
  period <- subset(period,id!="SH1_2013__CIFORcameroon")
  period <- UpdateSubsets(period,Dir=savedir) 
  # make an optional subset
  period2 <- subset(period,id=="SH1_2013__CIFORcameroon")
  
  # To download the pixels  
  k = 1
  MODISSubsets(LoadDat = period2[k:nrow(period2),],Products = product,Bands = bands,Size = pixel,SaveDir = savedir,StartDate = T)
  for(i in 1:length(unique(period$id))){
    period2 <- subset(period,id%in%(unique(period$id)[i]))
    print(paste("Nr.",i,"-",unique(period2$id)))
    try(a <- MODISSubsets(LoadDat = period2,Products = product,Bands = bands,Size = pixel,SaveDir = savedir,StartDate = T)
      )
    print(paste(length(unique(period$id))-i,"remaining..."))
  }
  
  
  #source("MODIS_YieldCalculation.R")
  #summaries <- MODIS_yield(savedir,period,c(0))
  
  #summaries$missing
  pred_gisdata <- readRDS("Data/PredictData.rds")
  i <- 1
  period <- data.frame(lat=coord[i:nrow(coord),1],long=coord[i:nrow(coord),2],start.date=dates[i:nrow(coord),1],end.date=dates[i:nrow(coord),2],id=sites.pred$SSS[i:nrow(sites.pred)])

  MODISSummaries(LoadDat = period, Product = "MOD13Q1",Dir = "PredictsMODIS2years_PREDICTS/",
                 Bands = c("250m_16_days_EVI","250m_16_days_NDVI"),
                 ValidRange = c(-2000,10000), NoDataFill = -3000, ScaleFactor = 0.0001,
                 StartDate = TRUE, QualityScreen = TRUE, QualityThreshold = 0,
                 QualityBand = "250m_16_days_pixel_reliability",
                 Interpolate=T,Mean=T,Yield=T)

  file.remove("PredictsMODIS2years_PREDICTS/MODIS_Summary_MOD13Q1_2014-09-25_h15-m38-s34.csv")
  file.remove("PredictsMODIS2years_PREDICTS/MODIS_Data_MOD13Q1_2014-09-25_h15-m38-s34.csv")
  mr <- read.csv("PredictsMODIS2years_PREDICTS/MODIS_Summary_MOD13Q1_2014-09-25_h17-m16-s7.csv",header=T)
  md <- read.csv("PredictsMODIS2years_PREDICTS/MODIS_Data_MOD13Q1_2014-09-25_h17-m16-s7.csv",header=T)  
  stopifnot(md$id==period$id)
  which(!(period$id%in%md$id))
  m_ev <- subset(mr,data.band=="250m_16_days_EVI")
  m_ev$start.date <- year(m_ev$start.date)
  m_ev$end.date <- year(m_ev$end.date)
  m_ev$index <- paste(m_ev$lat,m_ev$long,m_ev$start.date,m_ev$end.date,sep="_")
  period$index <- paste(period$lat,period$long,period$start.date,period$end.date,sep="_")
  m <-join(period,m_ev,by="index",match = "first")  
  pred_gisdata$meanEVI <- m$mean.band
  pred_gisdata$yieldEVI <- m$band.yield

  m_ev <- subset(mr,data.band=="250m_16_days_NDVI")
  m_ev$start.date <- year(m_ev$start.date)
  m_ev$end.date <- year(m_ev$end.date)
  m_ev$index <- paste(m_ev$lat,m_ev$long,m_ev$start.date,m_ev$end.date,sep="_")
  period$index <- paste(period$lat,period$long,period$start.date,period$end.date,sep="_")
  m <-join(period,m_ev,by="index",match = "first")
  pred_gisdata$meanNDVI <- m$mean.band
  pred_gisdata$yieldNDVI <- m$band.yield
  
  pairs.panels(pred_gisdata[,5:9])

  # Check for consistency with Jorns function
  # Call new script
  addModis <- read.csv("MODIS pixel summaries.txt",sep="\t",header=T)
  names(addModis) <- paste0("upp.",names(addModis));names(addModis)[1] <- "ID"
  sub <- subset(o,data.band=="250m_16_days_NDVI")
  o <- join(sub,addModis,by="ID")
  plot(o$mean.band,o$upp.mean.ndvi)

  ## ------------------------- ##
  #### FC Hansen ####
  fc2 <- raster("/media/GISbox/Raster/Landcover/HANSEN_2013/Treecover2000/fc2000.vrt")
  sp <- SpatialPointsDataFrame(cbind(sites.pred$Longitude,sites.pred$Latitude),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),data=pred_gisdata)
  pred_gisdata$FC2000 <- raster::extract(fc2,sp)
  rm(fc2)
  
  # ---------------- #
  #### BIOCLIM baseline data #### 
  unlink("/media/GISbox/Temp/*");dir.create("/media/GISbox/Temp/")
  unzip("/media/GISbox/Raster/Climate/WorldClim_Africa_Kite/Baseline/tbio_wc30s.zip",exdir = "/media/GISbox/Temp/")
  unzip("/media/GISbox/Raster/Climate/WorldClim_Africa_Kite/Baseline/mbio_wc30s.zip",exdir = "/media/GISbox/Temp/")
  bc <- lapply(list.files(path=paste0("/media/GISbox/Temp/"),pattern="*.tif",full.names=T),raster)
  bioclim <-stack(bc);rm(bc)  
  bcd <- as.data.frame(raster::extract(bioclim,sp))
  pred_gisdata <- cbind(pred_gisdata,bcd)
  unlink("/media/GISbox/Temp/*")
  file.remove("/media/GISbox/Temp")

  # Save everything
  saveRDS(pred_gisdata,"Data/PredictData.rds")
  
}
predicts_data_load <- function(sites.pred){
  # Load Container
  data <- readRDS("Data/PredictData.rds")
  out <- cbind(sites.pred,data)
  print("Appended Data. Now correcting NDVI for seasonality.")  
  sp <- subset(out,select=c("Latitude","Longitude","yieldNDVI","Biome","PREDICTS.LU","PREDICTS.LUI"))
  sp <- SpatialPointsDataFrame(cbind(sp$Longitude,sites.pred$Latitude),proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),data=sp)
  out$yield.ndvi.corr <- correctNDVIyield(sp,F)
  #yield.ndvi.corr1 <- correctNDVIyield(sp,F)
  #yield.ndvi.corr2 <- correctNDVIyield(sp,T)
  #plot(yield.ndvi.corr1~yield.ndvi.corr2)
  #lines( par()$usr[1:2], par()$usr[3:4] )
  rm(sp)
  return(out)
}


correctNDVIyield <- function(sp,africlim=T){
  if(africlim==F){
    print("Sampling from 5arcmin bioclim")
    unzip("../GIS/Climate/bio_5m_bil.zip",exdir = "/tmp/bioclim")
    b4 <- raster("/tmp/bioclim/bio4.bil")
    sp$temp.seas <- extract(b4,sp)
    rm(b4)
    b15 <- raster("/tmp/bioclim/bio15.bil")
    sp$precip.seas <- extract(b15,sp)
    rm(b15)
    unlink("/tmp/bioclim/")
  } else {
    print("Sampling from 5arcseconds bioclim")
    unzip("../GIS/Climate/mbio_wc30s.zip",exdir = "/tmp/bioclim")
    unzip("../GIS/Climate/tbio_wc30s.zip",exdir = "/tmp/bioclim")
    b4 <- raster("/tmp/bioclim/bio1_wc30s.tif")
    sp$temp.seas <- (extract(b4,sp)*.1) # * 0.1
    rm(b4)
    b15 <- raster("/tmp/bioclim/bio15_wc30s.tif")
    sp$precip.seas <- extract(b15,sp)
    rm(b15)
    unlink("/tmp/bioclim/")
    
  }
  
  # Only use primary vegetation-minimal use sites to try to remove any effects of land use
  sp <- as.data.frame(sp)
  pristine.sites<-droplevels(sp[((sp$PREDICTS.LU=="Primary Vegetation") & 
                                   (sp$PREDICTS.LUI=="Minimal use") & 
                                   (!is.na(sp$yieldNDVI))),])
  
  pristine.sites <- na.omit(pristine.sites)
  # Now compare combined models
  # (constrain to additive model, because interactive model behaves strangely,
  # probably too few data points)
  m0<-lm(log(yieldNDVI)~1,data=pristine.sites)
  m1<-lm(log(yieldNDVI)~poly(temp.seas,2),data=pristine.sites)
  m2<-lm(log(yieldNDVI)~poly(precip.seas,3),data=pristine.sites)
  m3<-lm(log(yieldNDVI)~poly(temp.seas,2)+poly(precip.seas,3),data=pristine.sites)
  
  # Compare model fits
  a = AIC(m0,m1,m2,m3)
  best_model <- get(rownames(a[which(a$AIC==min(a$AIC)),]))
  
  ####Apply correction to diversity data####
  yield.corr<-sp$yieldNDVI-exp(predict(best_model,newdata=sp,na.action="na.exclude"))
  return(yield.corr)
}


ModisTileDownload <- function(x){
  # Direct Download of MODIS tiles
  mtiles <- readOGR("/media/GISbox/Raster/VegetationIndex/modis_sinusoidal_shape/","modis_sinusoidal_grid_world")
  dp = "ftp://e4ftl01.cr.usgs.gov/MOLT/MOD13Q1.005/"
  ddir <- "/media/GISbox/Raster/VegetationIndex/MODIS_MOD13Q1_PREDICTS/"
  sin_coord <- spTransform(sp,CRSobj = CRS(proj4string(mtiles)))
  ov <- over(sin_coord,mtiles);ov$cat <- NULL
  ov$h <- as.character(ov$h)
  l <- str_length(as.character(ov$v))
  ov$v[which(l<2)] <- paste0("0",ov$v[which(l<2)])
  ov$start.date <- pred_gisdata$start.date;ov$start.month <- pred_gisdata$start.month
  ov$end.date <- pred_gisdata$end.date;ov$end.month <- pred_gisdata$end.month
  ov$Sample_start_earliest <-  pred_gisdata$Sample_start_earliest
  ov$Sample_end_latest <- pred_gisdata$Sample_end_latest
  uov <- unique(ov) # Calculate only for unique rows <- reduces dataset to 30
  sink(paste0("/media/GISbox/Raster/VegetationIndex/","modis_dl.sh"))
  for(i in 1:nrow(uov)) {
    tile = paste0("h",uov$h[i],"v",uov$v[i])
    s = paste0(uov$start.date[i],"-",uov$start.month[i],"-01")
    e = paste0(uov$end.date[i],"-",uov$end.month[i],"-",days_in_month(as.numeric(uov$end.month[i])))
    cat(paste("echo","Processing",tile,"from",s,"to",e))
    cat("\n")
    p = paste0("h",uov$h[i],"v",uov$v[i]) 
    cat("mkdir -p",(paste0(ddir,p)))    
    cat("\n")    
    f <- paste0("modis_download.py -r -p MOD13Q1.005 -t ",tile," -e ",s," -f ",e," MODIS_MOD13Q1_PREDICTS/",p)     
    cat(f)
    cat("\n")
    cat(paste("echo fininished downloading.",nrow(unique(uov))-i,"remaining"))
    cat("\n")
    rm(tile,f,p)
  }
  sink()  
  # Stop and download
  stop("Run Modis_download.py script and update files on GISBOX")
 
  # Continue  
  pred_gisdata$meanNDVI <- NA
  pred_gisdata$meanEVI <- NA  
  pred_gisdata$meanInterpNDVI <- NA
  pred_gisdata$meanInterpEVI <- NA
  pred_gisdata$yieldNDVI <- NA
  pred_gisdata$yieldEVI <- NA
  #### Mass NDVI Yield Sampler ####
  require(gdalUtils)
  uov <- unique(ov)
  for(i in 1:nrow(uov)){
    rn <- as.numeric(rownames(match_df(x = ov,y =  uov[i,]))) # Gives us the row names
    tile = paste0("h",uov$h[i],"v",uov$v[i])
    d = paste0(ddir,tile,"/")
    std <- paste0(uov$start.date[i],formatC(yday(uov$Sample_start_earliest[i]),width=3,flag="0")) 
    end <- paste0(uov$end.date[i],formatC(yday(uov$Sample_end_latest[i]),width=3,flag="0"))
    h <- grep(pattern = "xml",list.files(d,pattern = "*.hdf",full.names = T,ignore.case = F),value = T,invert = T)
    # Using the MODIS Orgtime and prestack function to subset the datasets needed
    o <- orgTime(h,nDays = 16,begin = std,end = end,pillow=15)
    vi <- preStack(files = h,timeInfo = o)
    cat("In Total",length(vi),"different modis frames")
    cat("\n")
    # Load in the raster layers
    ras_list <- list()
    #type = NDVI & EVI = 16bit, reliability = 8bit integer
    for(lay in vi){
      
      mod <- getSds(lay,method="gdal")
      ndvi <- raster(mod$SDS4gdal[1])#raster(readGDAL(mod$SDS4gdal[1]))
      evi <- raster(mod$SDS4gdal[2])#evi2 <- raster(readGDAL(mod$SDS4gdal[2]))
      rel <- raster(mod$SDS4gdal[12])#raster(readGDAL(mod$SDS4gdal[12]))
      ras_list[lay] <- stack(ndvi,evi,rel)
      rm(ndvi,evi,rel,mod)
    }    
    # Seperate loop for individual coordinates
    for(p in rn){
      # Spatial subset and transformation
      s <- spTransform(sp[p,],CRS(proj4string(mtiles)))
      res = data.frame(NDVI=numeric(),EVI=numeric(),REL=numeric(),stringsAsFactors = F) 
      for(lay in vi){
        ex <- as.data.frame(raster::extract(ras_list[[lay]],s));names(ex) <- c("NDVI","EVI","REL")
        res <- rbind(res,ex)
        rm(ex)
      }
      names(res) <- c("NDVI","EVI","REL")
      res$NDVI <- res$NDVI *0.0001 # Apply Scaling Facor
      res$EVI <- res$EVI *0.0001 # Apply Scaling Facor
      # Kick out values with -1,23 
      cat("Data not useable from ",length(which(res$REL==0))," Dates out of",nrow(res))
      cat("\n")
      dates <- o$inputLayerDates#[which(res$REL==1)]
      #res <- res[which(res$REL==1),]
      
      # Calculate mean of EVI and NDVI
      pred_gisdata$meanNDVI[p] <-  mean(res$NDVI*0.0001)
      pred_gisdata$meanEVI[p] <-  mean(res$EVI*0.0001)
      if(nrow(res)>1&length(which(!is.na(res$NDVI)))>1){
        ######################
        # Linear interpolation and yield calculation
        ######################
        sout = approx(x=dates, y=res$NDVI*0.0001, method = "linear", n = ((max(dates[!is.na(res$NDVI)])- min(dates[!is.na(res$NDVI)]))-1))
        minobsndvi = min(res$NDVI*0.0001, na.rm = TRUE)    # minimum NDVI observed
        maxobsndvi = max(res$NDVI*0.0001, na.rm = TRUE)    # maximum NDVI observed
        ndvi.yield = (sum(sout$y) - minobsndvi*length(sout$x)) / length(sout$x) #(((365*length(years))-16)*365) # average annual yield  (i.e. work out daily yield * 365
        pred_gisdata$yieldNDVI[p] <- ndvi.yield
        pred_gisdata$meanInterpNDVI[p] <-mean(sout$y)
        
        minobsevi = min(res$EVI*0.0001, na.rm = TRUE)    # minimum NDVI observed
        maxobsevi = max(res$EVI*0.0001, na.rm = TRUE)    # maximum NDVI observed
        evout = approx(x=dates, y=res$EVI*0.0001, method = "linear", n = ((max(dates[!is.na(res$EVI)])- min(dates[!is.na(res$EVI)]))-1))
        evi.yield = (sum(evout$y) - minobsevi*length(evout$x)) / length(evout$x) #(((365*length(years))-16)*365) # average annual yield  (i.e. work out daily yield * 365
        pred_gisdata$yieldEVI[p] <- evi.yield
        pred_gisdata$meanInterpEVI[p] <- mean(evout$y)
        cat("Interpolated and extracted Data. yieldEVI =",evi.yield)
        cat("\n")        
      }
    }
    cat("Finished tile",tile)
    cat("\n")
    rm(ras_list)
  }
  
  cat("NA in",length(which(is.na(pred_gisdata$meanNDVI))),"out of",nrow(pred_gisdata),"\n")
  a = difftime(sites.pred$Sample_end_latest,sites.pred$Sample_start_earliest,units = "days")
  plot(as.vector(a)~pred_gisdata$yieldEVI)
  summary(l<-lm(as.vector(a)~pred_gisdata$yieldEVI));abline(l,col="red",lwd=2)
  setwd(ori.dir)
  
  
}

AppendHYDE <- function(Long,Lat,Year) {
  hyde <- raster("../GIS/HYDE/yoc30cropgrasuoppALLYRS.asc")
  proj4string(hyde) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  # Replace all 9999
  hyde[hyde == 9999] <- NA
  # Attached is the HYDE-derived estimate of the year when cells 
  # first became 30% converted to human land uses (cropland, pasture, urban).
  # Until 1000AD, the data are analyzed at 1000-year intervals, 
  # then 100-year intervals until 1900, 10-year intervals between 1900 and 2000,
  # with a final estimate in 2005. 
  # A value of 9999 means that a cell is estimated never to have been under 30% human land use.  
  
  # turn this data frame into a Spatial Point object
  coords <- as.data.frame(cbind(Long, Lat))  
  coordinates(coords) <- ~Long + Lat
  proj4string(coords) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  
  # Extract from HYDE
  res <- extract(hyde,coords)
  res <- (Year - res)
  res[which(res < 0)] <- 0
  #res[which(is.na(res))] <- 0
  
  return(res)  
}


HYDEcustomFigure <- function(sites,sites.pred){
  hyde <- raster("../GIS/HYDE/yoc30cropgrasuoppALLYRS.asc")
  proj4string(hyde) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  # Replace all 9999
  hyde[hyde == 9999] <- NA
  # Attached is the HYDE-derived estimate of the year when cells 
  # first became 30% converted to human land uses (cropland, pasture, urban).
  # Until 1000AD, the data are analyzed at 1000-year intervals, 
  # then 100-year intervals until 1900, 10-year intervals between 1900 and 2000,
  # with a final estimate in 2005. 
  # A value of 9999 means that a cell is estimated never to have been under 30% human land use.  
  
  # First PREDICTS wide
  sites.pred$Source_ID <- droplevels(sites.pred$Source_ID)
  sites.pred$SS <- droplevels(sites.pred$SS)
  sp1 <- data.frame(id=(sites.pred$SS),Long=sites.pred$Longitude,Lat=sites.pred$Latitude)
  coordinates(sp1) <- ~Long + Lat
  proj4string(sp1) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  
  ctrs <- lapply( unique(sp1$id), function(x) gCentroid(sp1[which(sp1$id==x),],byid = F ) )
  ctrs <- setNames( ctrs , unique(sp1$id ) )
  res <- data.frame(SS=names(ctrs),ex=NA,Year=NA,Est=NA)
  res$Year <- unlist(lapply( unique(res$SS), function(x) min(unique(year(sites.pred$Sample_start_earliest[which(sites.pred$SS==x)]))) ) )

  for(SS in names(ctrs)){
    s <- ctrs[[SS]]
    res$ex[which(res$SS==SS)] <- extract(hyde,s)
    res$Est <- (res$Year - res$ex)
    # res$Est[which(res$Est < 0)] <- 0 # Set sites that haven't been altered before the beginning of sampling to 0    
  }
  
  # And field data
  s <- subset(sites,Grouping=="Aves")
  sp1 <- data.frame(id=(s$Transect),Long=s$Long,Lat=s$Lat)
  coordinates(sp1) <- ~Long + Lat
  proj4string(sp1) <- CRS("+proj=longlat +datum=WGS84 +no_defs")  
  ctrs <- lapply( unique(sp1$id), function(x) gCentroid(sp1[which(sp1$id==x),],byid = F ) )
  ctrs <- setNames( ctrs , unique(sp1$id ) )
  res2 <- data.frame(SS=names(ctrs),ex=NA,Year=2014,Est=NA)
  
  for(SS in names(ctrs)){
    s <- ctrs[[SS]]
    res2$ex[which(res2$SS==SS)] <- extract(hyde,s)
    res2$Est <- (res2$Year - res2$ex)
    # res$Est[which(res$Est < 0)] <- 0 # Set sites that haven't been altered before the beginning of sampling to 0    
  }
  res$Type <- "PREDICTS"
  res2$Type <- "Fieldata"
  r <- rbind(res,res2)
  b <- ddply(predictsdata,("SS"),summarise,
             SpeciesRichness = length(unique(Best_guess_binomial))
  )
  r <- join(r,b)
  library(ggthemr)
  ggthemr("grape",type = "outer")
  g <- qplot(Est,fill=Type,data=r[-19,],binwidth=3,xlab="Time since 30% anthropogenic conversion (years)",main = paste("Based on ",length(r$SS[-19])," studies (SS)"))
  ggsave("~/HYDE_estimated_WithoutLachat.png",g,dpi = 400)
  
  ddply(r,.(Type),summarise,Avg=mean(Est,na.rm = T),SD=sd(Est,na.rm=T))
  
}
