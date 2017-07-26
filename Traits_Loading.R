#### Package and Data loading ####
#rm(list=ls())
library(marfunky)
standardPackages(c("Standard","Plot","Analysis"))
library(yarg)
library(roquefort)
library(letsR)
library(taxize)
library(rgbif)
library(spocc)
setwd("~/Documents/Studium/Kopenhagen/Masters/Analysis_PREDICTS")
# Load scientific names of my field data
species <- read.xls("../FieldData/Species.xls",sheet = 1)[,3]
# And Dickens data
speciesD <- read.xls("../FieldData/Dickens_Taita_Kilimanjaro_Stock_Density.xls",sheet=1)
rownames(speciesD) <- speciesD$PlotNumb
speciesD <- speciesD[,c(-1,-ncol(speciesD))] # Drop last and first column

# Predicts data output from 28-07-2014
predictsdata <- readRDS("Data/martin-2014-07-28-16-47-26.rds")

predictsdata <- DropInvalidMetricsAndMethods(predictsdata)
predictsdata <- CorrectSamplingEffort(predictsdata)

ind = which(!(predictsdata$Species)=="")
ind = 1:nrow(predictsdata)
pd_spec <- predictsdata$Best_guess_binomial[ind]
#pd_spec_group <- predictsdata$Higher_taxon[ind]
#rm(predictsdata,ind)
spec_lookup <- unique(pd_spec[which(!(pd_spec==""))])

pd <- predictsdata %>% filter(Higher_taxon=="Aves")
(unique(pd$Best_guess_binomial)) #564
length(which(pd$Best_guess_binomial!="")) / nrow(pd)

#### IUCN Data ####
# My data
iucn_sum <- letsR::lets.iucn(species,F) # done
iucn_ha <- letsR::lets.iucn.ha(species,F)
iucn_hi <- letsR::lets.iucn.his(species,F)

iucn_field <- join(iucn_sum,iucn_hi,by=,match = "all")
iucn_field <- cbind(iucn_field,iucn_ha)
write.csv(iucn_field,"TraitsData/IUCN_field_data.csv",row.names=F)


# Alternative IUCN query!
# Match with Full data
source("isForestSpecialist.R")
data <- read.csv("IUCN_SpeciesMaster/all.csv")
d <- data.frame(species=as.character(species))
d$speciesID <- data$Red.List.Species.ID[match(d$species,data$Scientific.Name)]
d$ForestSpec<- sapply(d$speciesID,isForestSpecialist)
saveRDS(d,"newIUCNHabitatinfo.rds")
rm(d)

# Dickens data
sdn <- names(speciesD)
sdn <- str_replace(sdn,pattern = "_",replacement = " ")

iucn_sum <- letsR::lets.iucn(sdn,F) # done
iucn_ha <- letsR::lets.iucn.ha(sdn,F)
iucn_hi <- letsR::lets.iucn.his(sdn,F)
iucn_field <- join(iucn_sum,iucn_hi,by=,match = "all")
iucn_field <- cbind(iucn_field,iucn_ha)
write.csv(iucn_field,"TraitsData/IUCN_Dfield_data.csv",row.names=F)

a <- read.csv("TraitsData/IUCN_field_data.csv")
head(a)
rm(a)
# PREDICTS data
iucn_sum1 <- letsR::lets.iucn(spec_lookup[1:1500],count = F)
write.csv(iucn_sum1,"TraitsData/IUCN_PREDICTS1.csv",row.names=F)
iucn_sum2 <- letsR::lets.iucn(spec_lookup[1500:2500],F)
write.csv(iucn_sum2,"TraitsData/IUCN_PREDICTS2.csv",row.names=F)
iucn_sum3 <- letsR::lets.iucn(spec_lookup[2500:length(spec_lookup)],F)
write.csv(iucn_sum3,"TraitsData/IUCN_PREDICTS3.csv",row.names=F)
# Habitat data
iucn_ha1 <- letsR::lets.iucn.ha(spec_lookup[1:1000],count = F)
write.csv(iucn_ha1,"TraitsData/IUCN_PREDICTS_ha1.csv",row.names=F)
iucn_ha2 <- letsR::lets.iucn.ha(spec_lookup[1000:1500],F)
write.csv(iucn_ha2,"TraitsData/IUCN_PREDICTS_ha2.csv",row.names=F)
iucn_ha3 <- letsR::lets.iucn.ha(spec_lookup[1500:2000],F)
write.csv(iucn_ha3,"TraitsData/IUCN_PREDICTS_ha3.csv",row.names=F)
iucn_ha4 <- letsR::lets.iucn.ha(spec_lookup[2000:2500],F)
write.csv(iucn_ha4,"TraitsData/IUCN_PREDICTS_ha4.csv",row.names=F)
iucn_ha5 <- letsR::lets.iucn.ha(spec_lookup[2500:3000],F)
write.csv(iucn_ha5,"TraitsData/IUCN_PREDICTS_ha5.csv",row.names=F)
iucn_ha6 <- letsR::lets.iucn.ha(spec_lookup[3000:length(spec_lookup)],F)
write.csv(iucn_ha6,"TraitsData/IUCN_PREDICTS_ha6.csv",row.names=F)


# Alternative approach
source("isForestSpecialist.R")
data <- read.csv("IUCN_SpeciesMaster/all.csv")
d <- data.frame(species=as.character(spec_lookup))
d$speciesID <- data$Red.List.Species.ID[match(d$species,data$Scientific.Name)]
d$ForestSpec<- sapply(d$speciesID,isForestSpecialist)
saveRDS(d,"newIUCNHabitatinfo_pred.rds")
rm(d)


# Load them all and save as rds
h1 <- read.csv("TraitsData/IUCN_PREDICTS_ha1.csv",header=T)
h2 <- read.csv("TraitsData/IUCN_PREDICTS_ha2.csv",header=T)
h3 <- read.csv("TraitsData/IUCN_PREDICTS_ha3.csv",header=T)
h4 <- read.csv("TraitsData/IUCN_PREDICTS_ha4.csv",header=T)
h5 <- read.csv("TraitsData/IUCN_PREDICTS_ha5.csv",header=T)
h6 <- read.csv("TraitsData/IUCN_PREDICTS_ha6.csv",header=T)
h <- rbind(h1,h2,h3,h4,h5,h6);rm(h1,h2,h3,h4,h5,h6)
h <- h[which(duplicated(h)==F),]
d <- data.frame(Species=spec_lookup,ha=h)
saveRDS(h,"TraitsData/IUCN_PREDICTS_Habitats.rds")

#### Make a tile of percentual contribution of each data ####
pred_threat <- read.csv("TraitsData/IUCN_PREDICTS_Threats.csv",header=T)
pred_threat$Best_guess_binomial <- pred_threat$Species
pred_threat <- unique(pred_threat)
names(pred_threat)[3] <- "IUCN_Status"
table(pred_threat$IUCN_Status)
pred_threat$IUCN_cont <- lets.iucncont(pred_threat$IUCN_Status)

pred_threat$Species[which(!pred_threat$Species%in%spec_lookup)]

predictsdata2 <- merge(predictsdata,pred_threat,by = "Best_guess_binomial",all.x = T)
# Subset only to Birds
d <- data.frame(from=unique(predictsdata2$Study_common_taxon),
                to=c("Invertebrates","Invertebrates","Invertebrates","Birds","Reptiles","Mammals","Other",
                     "Invertebrates","Plants","Amphibia","Plants","Plants","Plants","Plants","Mammals","Other","Plants","Amphibia",
                     "Mammals","Invertebrates")) # for matching
predictsdata2$Grouping <- NA
predictsdata2$Grouping <-  d$to[match(predictsdata2$Study_common_taxon,d$from)]
rm(d)
predictsdata2 <- subset(predictsdata2,Grouping=="Plants")

# Kick out unsuitable studies
pred <- subset(predictsdata2,Source_ID!="GP1_2012__Strauch")
pred <- subset(pred,!(Predominant_habitat%in%c("Pasture")))
pred <- subset(pred,year(pred$Sample_start_earliest)>=2000)

# Reclassify LandUse values
pred$LandUse<-paste(pred$Predominant_habitat)
pred$LandUse[which(pred$LandUse=="Primary forest")]<-"Primary Vegetation"
pred$LandUse[which(pred$LandUse=="Primary non-forest")]<-"Primary Vegetation"
pred$LandUse[which(pred$LandUse=="Young secondary vegetation")]<-"Secondary Vegetation"
pred$LandUse[which(pred$LandUse=="Intermediate secondary vegetation")]<-"Secondary Vegetation"
pred$LandUse[which(pred$LandUse=="Mature secondary vegetation")]<-"Secondary Vegetation"
pred$LandUse[which(pred$LandUse=="Secondary vegetation (indeterminate age)")]<-"Secondary Vegetation"
pred$LandUse[which(pred$LandUse=="Secondary non-forest")]<-"Secondary Vegetation"
pred$LandUse[which(pred$LandUse=="Cannot decide")]<-NA
pred$LandUse<-factor(pred$LandUse)
pred$LandUse<-relevel(pred$LandUse,ref="Primary Vegetation")
pred$PREDICTS.LU <- pred$LandUse
pred$PREDICTS.LUI <- pred$Use_intensity


(mat <- prop.table(table(pred$IUCN_Status,pred$PREDICTS.LU),2))
# Make a Tile plot

m <- melt(mat,na.rm = T)
m$Var2 <- ordered(m$Var2,c("Primary Vegetation", "Secondary Vegetation","Plantation forest","Pasture", "Cropland","Urban") )
m$Var1 <- ordered(m$Var1,c("EX","EW","CR","EN","VU","NT","LC","DD","NE"))

theme_set(theme_classic(base_size=12,base_family = "sans"))
g <- ggplot(m, aes(Var1, Var2, fill = value))
g <- g + geom_tile() + coord_flip()
g <- g + geom_text(aes(fill = m$value, label = paste0(round(m$value*100, 2),"%")),col="white")
g <- g + scale_fill_continuous(guide="none") 
g <- g + ggtitle(" PREDICTS IUCN current Status (2014)") + xlab("") + ylab("")
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g
ggsave("~/PREDICTS_IUCN.png",plot=g,dpi=400)

site_threat <- read.csv("TraitsData/IUCN_field_data.csv",header=T)

source("../FieldData/PrepareFiles.R")
sp <- species
rownames(sp) <- sites$PREDICTS.LU[which(sites$SiteName%in%rownames(sp))]
colnames(sp) <- site_threat$Status
m <- melt(sp)

ag <- aggregate(list(value=m$value),by=list(Var1=m$Var1,Var2=m$Var2),FUN=sum)
m <- ddply(ag, .(Var1), summarise, Var2 = Var2, value = value / sum(value))
m$Var1 <- ordered(m$Var1,c("Primary Vegetation", "Secondary Vegetation","Plantation forest","Pasture", "Cropland","Urban") )
m$Var2 <- ordered(m$Var2,c("EX","EW","CR","EN","VU","NT","LC","DD","NE"))
theme_set(theme_classic(base_size=12,base_family = "sans"))
g <- ggplot(m, aes(Var1, Var2, fill = value))
g <- g + geom_tile()
g <- g + geom_text(aes(fill = m$value, label = paste0(round(m$value*100, 2),"%")),col="white")
g <- g + scale_fill_continuous(guide="none") 
g <- g + ggtitle(" Fieldwork IUCN current Status (2014)") + xlab("") + ylab("")
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g
ggsave("~/Fieldwork_IUCN.png",plot=g,dpi=400)


# Numerical comparison - Change to Plant above 
site_threat <- read.csv("TraitsData/IUCN_field_data.csv",header=T)
site_threat <- read.csv("TraitsData/IUCN_Dfield_data.csv",header=T)

ag1 <- aggregate(list(pred=pred$IUCN_cont),by=list(pred$PREDICTS.LU),FUN = function(x)mean(x,na.rm=T))
#ag1 <- rbind(ag1,data.frame(Group.1="Urban",pred=NA))
sub <- subset(sites,Grouping=="Plantae")
sp <- as.matrix(plantdata)#species
rownames(sp) <- sub$PREDICTS.LU[which(sub$SiteName%in%rownames(sp))]
colnames(sp) <- lets.iucncont(site_threat$Status)
m <- melt(sp)
m <- m[which(!is.na(m$Var2)),]

ag2 <- aggregate(list(field=m$value),by=list(m$Var1),FUN=function(x)mean(x,na.rm=T))
agg <- join(ag1,ag2)

t.test(agg$pred,agg$field,exact = T,paired = T)


# Change Threat category over time since 2000
# Community-wide mean threat rank
# Does Threat rank correspond to environmental factors / land-use intensity

#### Range of occurence ####
# Read in GBIF 1 Degree cell id shapefile
library(rgdal)
grid <- readOGR("../GIS/","GBIF_Africa1degCells")

# d <- occ(as.character(species[1]),limit = 1000000,from = "gbif",gbifopts = list(hasCoordinate = TRUE,continent = "africa"))
# sp <- occ2sp(d)
# out <- density_spplist(taxonconceptKey = key,spplist = "none",originisocountrycode = c("KE","TZ"),listcount='counts')
# dismo::gbif(genus ="Melaenornis", species = "fischeri")
# splist <- as.character(species[101:103])
# out <- occurrencelist_many(splist,originisocountrycode =c("KE","TZ")  ,coordinatestatus = TRUE,maxresults = 1000000)
# gbifmap_list(out,region =c("Tanzania","Kenya") )

# For my species
res <- data.frame()
for(id in species){
  print(paste("Processing ",id))
  key <- name_suggest(q=as.character(id), rank='species')$key[1] #Best guess
  o <- densitylist(taxonconceptKey = key) # get number of occupied cells
  o <- o[which(o$cellid%in%grid$cellid),] # subset to those in Africa
  if(length(o$cellid)==0)print("Species not found (in Africa)")
  d <- data.frame(SpeciesID=id,NrCells=length(o$cellid),OccupiedCells=I(list(o$cellid)))
  res <- rbind(res,d)
  rm(o,d)
  print("------")
}
saveRDS(res,"TraitsData/GBIF_OCC-Fielddata.rds")

# For Dickens species
resD <- data.frame()
for(id in sdn){
  print(paste("Processing ",id))
  key <- name_suggest(q=as.character(id), rank='species')$key[1] #Best guess
  o <- densitylist(taxonconceptKey = key) # get number of occupied cells
  o <- o[which(o$cellid%in%grid$cellid),] # subset to those in Africa
  if(length(o$cellid)==0)print("Species not found (in Africa)")
  d <- data.frame(SpeciesID=id,NrCells=length(o$cellid),OccupiedCells=I(list(o$cellid)))
  resD <- rbind(resD,d)
  rm(o,d)
  print("------")
}
saveRDS(resD,"TraitsData/GBIF_OCC-DFielddata.rds")

# For PREDICTS
spec_lookup <- unique(pd_spec[which(!(pd_spec==""))])
res <- data.frame()
for(id in spec_lookup){
  print(paste("Processing ",id))
  key <- name_suggest(q=as.character(id), rank='species')$key[1] #Best guess
  #o <- occ_search(taxonKey=key, limit=20,hasCoordinate = T,geometry = bbox2wkt(bbox = bbox(grid)), return='data')
  ## Deprecated
  o <- densitylist(taxonconceptKey = key) # get number of occupied cells
  o <- o[which(o$cellid%in%grid$cellid),] # subset to those in Africa
  if(length(o$cellid)==0)print("Species not found (in Africa)")
  d <- data.frame(SpeciesID=id,NrCells=length(o$cellid),OccupiedCells=I(list(o$cellid)))
  #gbifmap_dens(o)
  res <- rbind(res,d)
  rm(o,d)
  print(paste(length(spec_lookup)-nrow(res), "names remaining."))
  print("------")
}
saveRDS(res,"TraitsData/GBIF_OCC-PREDICTS.rds")


# Plot a bit
occ_pd <- readRDS("TraitsData/GBIF_OCC-PREDICTS.rds")
occ_pd$NrCells[which(occ_pd$NrCells==0)] <- NA
hist(occ_pd$NrCells,n=50,col="grey")
abline(v=median(occ_pd$NrCells,na.rm = T),col="blue",lwd=2)

# wide-ranged and narrow-ranged
# species whose area of occupancy exceeded the median for the broad taxonomic group
# were classed as wide-ranged and the others as narrow-ranged
d <- data.frame(from=unique(predictsdata$Study_common_taxon),
                to=c("Invertebrates","Invertebrates","Invertebrates","Birds","Reptiles","Mammals","Other",
                     "Invertebrates","Plants","Amphibia","Plants","Plants","Plants","Plants","Mammals","Other","Plants","Amphibia",
                     "Mammals","Invertebrates")) # for matching
predictsdata$Grouping <- NA
predictsdata$Grouping <-  d$to[match(predictsdata$Study_common_taxon,d$from)]
rm(d)

predictsdata$NrCells <- occ_pd$NrCells[match(predictsdata$Best_guess_binomial,occ_pd$SpeciesID)]
predictsdata$DistBound <- occ_pd$DistBound[match(predictsdata$Best_guess_binomial,occ_pd$SpeciesID)]

OvMed <- ddply(predictsdata,.(Grouping),summarise,
      med = median(NrCells,na.rm = T),
      lowq = quantile(NrCells,probs = .25,na.rm = T),
      highq = quantile(NrCells,probs = .75,na.rm = T))
# Fucking big one-line
predictsdata$Range <- ifelse(predictsdata$NrCells>OvMed$med[match(predictsdata$Grouping,OvMed$Grouping)],"wide-ranged","narrow-ranged")
# The alternative way
occ_pd$DistBound <- sapply(occ_pd$OccupiedCells,returnExtentPerimeter)  

hist(occ_pd$DistBound,n=100)
abline(v=median(occ_pd$DistBound,na.rm = T),col="blue",lwd=2)

occ_fw <- readRDS("TraitsData/GBIF_OCC-Fielddata.rds")
occ_fw$NrCells[which(occ_fw$NrCells==0)] <- NA
occ_fw$Range <- ifelse(occ_fw$NrCells>median(occ_fw$NrCells,na.rm = T),"wide-ranged","narrow-ranged")

hist(occ_fw$NrCells,n=50,col="grey")
abline(v=median(occ_fw$NrCells,na.rm = T),col="blue",lwd=2)
# Alternative of classifying the median. Calculate great circle distance between bbox corners
returnExtentPerimeter <- function(cid){ 
  cid <- unlist(cid)
  if(length(cid)==0){
    return(NA)
  } else {
    print(paste("Processing ",length(cid),"tiles"))  
    sg <- subset(grid,cid%in%cellid)
    b <- bbox(sg)
    d <- sp::spDists(t(b),longlat = T)[1,2]
    print(class(d))
    return(as.numeric(d))
  }
}

occ_fw$DistBound <- sapply(occ_fw$OccupiedCells,returnExtentPerimeter)                          
hist(occ_fw$DistBound,n=100)
abline(v=median(occ_fw$DistBound,na.rm = T),col="blue",lwd=2)
occ_fw$Range2 <- ifelse(occ_fw$DistBound>median(occ_fw$DistBound,na.rm = T),"wide-ranged","narrow-ranged")

plot(predictsdata$Grouping,predictsdata$NrCells)


#### Prepare Traits per site ####
source("FieldData/PrepareFiles.R")

repL <- function(mat){
  m2 <- mat
  for (r in seq(nrow(m2))){
    for (c in seq(ncol(m2))){
      if(m2[r,c]>0) m2[r, c] <- colnames(m2)[c] else m2[r, c] <- NA
    }
  }
  return(m2)
}

# ----------------------- #
# First for all field sites including Dickens
print("Percentage of narrow-ranged species")

occ_f <- readRDS("TraitsData/GBIF_OCC-Fielddata.rds")
occ_f2 <- readRDS("TraitsData/GBIF_OCC-DFielddata.rds")
occ <- rbind(occ_f,occ_f2);rm(occ_f,occ_f2)
occ$NrCells[which(occ$NrCells==0)] <- NA
occ$Range <- ifelse(occ$NrCells>median(occ$NrCells,na.rm = T),"wide-ranged","narrow-ranged")
getRange <- function(x){
  if(is.na(x))return(NA)else
    x <- str_replace(x,"_"," ")
    y = occ$Range[which(occ$SpeciesID==x)]
  if(is.na(y))return("Unknown")else return(y)  
}
spl <- read.xls("FieldData/Species.xls",sheet = 1,perl="C:/Perl64/bin/perl.exe")
colnames(species) <- spl$Scientific.Name[which(colnames(species)%in%spl$Common.Name)];rm(spl)
sp_c <- species;sp_c <- repL(sp_c)
sp_c <- apply(sp_c,1:2,getRange)
val <- apply(sp_c,1,FUN=function(x) length(which(x=="narrow-ranged")))/ apply(sp_c,1,function(x) length(which(!is.na(x))))

sp_c <- speciesD;sp_c <- repL(sp_c)
sp_c <- apply(sp_c,1:2,getRange)
val2 <- apply(sp_c,1,FUN=function(x) length(which(x=="narrow-ranged")))/ apply(sp_c,1,function(x) length(which(!is.na(x))))

sites$SpecINF_NarrowRanged <- c(val,val2)

# Alternative
# Calculate Community-wide mean range size
sp_c <- species; spc_c <- repL(sp_c)
sp_c[sp_c==0] <- NA
spc2 <- melt(sp_c,na.rm = T)
spc2 <- merge(spc2,occ,by.x ="Var2",by.y="SpeciesID")
d1 <- ddply(spc2,.(Var1),summarise,AvgRange=mean(log(NrCells),na.rm=T))

sites$SpecINF_AvgRange <- d1$AvgRange[match(sites$SiteName,d1$Var1)]

rm(sp_c,spc2,d1,val,val2,getRange)

# Last alternative
# Get traits from birdlife international
brs <- readRDS("BirdLifeRangeSizes.rds")
sp_c <- species; spc_c <- repL(sp_c)
sp_c[sp_c==0] <- NA
spc2 <- melt(sp_c,na.rm = T)
#library(taxize)
#d <- data.frame(spc2$Var2)
#spc2[which(is.na(match(spc2$Var2,brs$SCINAME))),]
spc2 <- merge(spc2,brs,by.x ="Var2",by.y="SCINAME")
d1 <- ddply(spc2,.(Var1),summarise,AvgRange=mean(log1p(DegOccupiedArea),na.rm=T))
sites$SpecINF_AvgRange2 <- d1$AvgRange[match(sites$SiteName,d1$Var1)]

rm(sp_c,spc2,d1,val,val2,getRange)

# ----- NEXT ---- #
print("Percentage of specialists ")
f1 <- read.csv("TraitsData/IUCN_field_data.csv")
f2 <- read.csv("TraitsData/IUCN_Dfield_data.csv")
field <- rbind(f1,f2);rm(f1,f2)

ha <- field[,c(44:51,57)];ha <- ifelse(ha==0,NA,1)
rs <- rowSums(ha,na.rm = T);rm(ha)
rs[rs==0] <- NA # If not occuring in at least one terrestrial habitat set to unknown
rs <- ifelse(rs==1,"habitat-specialist","habitat-generalist")
d <- data.frame(Species=as.character(field$Species),ha=as.character(rs));rm(rs)
getHa <- function(x) if(!is.na(x)) return(as.character(d$ha[which(d$Species==x)])) else return(NA)
sp_c <- species;sp_c <- repL(sp_c)
sp_c <- apply(sp_c,1:2,getHa)
val <- apply(sp_c,1,FUN=function(x) length(which(x=="habitat-specialist")))/ apply(sp_c,1,function(x) length(which(!is.na(x))))

sp_c <- speciesD
colnames(sp_c) <- str_replace(colnames(sp_c),"_"," ")
sp_c <- repL(sp_c)
sp_c <- apply(sp_c,1:2,getHa)
val2 <- apply(sp_c,1,FUN=function(x) length(which(x=="habitat-specialist")))/ apply(sp_c,1,function(x) length(which(!is.na(x))))

sites$SpecINF_HabitatSpec <- c(val,val2)
rm(sp_c,val,val2,getHa,d)

# ALTERNATIVE 
# Do the same with forest specialists!
f <- readRDS("newIUCNHabitatinfo.rds")
d <- data.frame(Species=as.character(f$species),ha=as.character(ifelse(f$ForestSpec=="Forest-specialist","Forest species","Non-forest species")))
getHa <- function(x) if(!is.na(x)) return(as.character(d$ha[which(d$Species==x)])) else return(NA)
sp_c <- species;sp_c <- repL(sp_c)
sp_c <- apply(sp_c,1:2,getHa)
val <- apply(sp_c,1,FUN=function(x) length(which(x=="Forest species")))/ apply(sp_c,1,function(x) length(which(!is.na(x))) )

sites$SpecINF_ForestSpec <- c(val,rep(NA,nrow(speciesD)))

# ----- NEXT ---- #
lets.iucncont2 <- function (x, dd = NA, ne = NA) {
  x <- as.matrix(x)
  for (i in 1:ncol(x)) {
    if (is.factor(x[, i])) {
      x[, i] <- as.numeric(levels(x[, i]))[x[, i]]
    }
  }
  x[(x == "EX" | x == "EW")] <- 5
  x[x == "CR"] <- 4
  x[x == "EN"] <- 3  
  x[x == "VU"] <- 2
  x[x == "NT"] <- 1
  x[x == "LC"] <- 0
  x[x == "DD"] <- dd
  x[x == "NE"] <- ne
  if (ncol(x) == 1) {
    x <- as.numeric(as.vector(x))
  }
  else {
    x <- as.data.frame(x)
  }
  return(x)
}
print("Average community-wide threat score")
tr <- lets.iucncont2(field[,3],dd = NA,ne = NA)
d <- data.frame(Species=as.character(field$Species),tr=tr);rm(tr)
getTr <- function(x) if(!is.na(x)) return(as.numeric(d$tr[which(d$Species==x)])) else return(NA)
sp_c <- species;sp_c <- repL(sp_c)
sp_c <- apply(sp_c,1:2,getTr)
val <- rowMeans(sp_c,na.rm = T)

sp_c <- speciesD
colnames(sp_c) <- str_replace(colnames(sp_c),"_"," ")
sp_c <- repL(sp_c)
sp_c <- apply(sp_c,1:2,getTr)
val2 <- rowMeans(sp_c,na.rm = T)

sites$SpecINF_IUCNStatus <- c(val,val2)
rm(sp_c,val,val2,getTr,d)

# Alternative:
# Proportion of protected species
tr <- lets.iucncont2(field[,3],dd = NA,ne = NA)
d <- data.frame(Species=as.character(field$Species),tr=tr)
getTr <- function(x) if(!is.na(x)) return(as.numeric(d$tr[which(d$Species==x)])) else return(-1)
sp_c <- species;sp_c <- repL(sp_c)
sp_c <- apply(sp_c,1:2,getTr)
val <- apply(sp_c,1,FUN=function(x) length(which(x >= 1))) / apply(sp_c,1,function(x) length(which(x >= 0)) )

sp_c <- speciesD
colnames(sp_c) <- str_replace(colnames(sp_c),"_"," ")
sp_c <- repL(sp_c)
sp_c <- apply(sp_c,1:2,getTr)
val2 <- apply(sp_c,1,FUN=function(x) length(which(x >= 1))) / apply(sp_c,1,function(x) length(which(x >= 0)))

sites <- subset(sites,Grouping=="Aves")
sites$SpecINF_IUCNStatus2 <- c(val)

cat("field sites finished!!!\n")

#----- now for PREDICTS#
print("Based on occurence")
occ <- readRDS("TraitsData/GBIF_OCC-PREDICTS.rds")
occ$NrCells[which(occ$NrCells==0)] <- NA
occ$Range <- ifelse(occ$NrCells>median(occ$NrCells,na.rm = T),"wide-ranged","narrow-ranged")
getRange <- function(x){
  if(x %in% occ$SpeciesID)
  { if(is.na(x))return(NA)else
    x <- str_replace(x,"_"," ")
  y = occ$Range[which(occ$SpeciesID==x)]
  #print(paste(x,"-",y))
  if( length(y)==0) return(NA) else{
    if(is.na(y)) return("Unknown") else return(y)
  }
  } else { return(NA)}
}
tr <- read.csv("TraitsData/IUCN_PREDICTS_Threats.csv",header=T)
tr$cont <- lets.iucncont2(tr$Status,dd = NA,ne = NA)
getTr <- function(x){ if(!is.na(x)) {
  if(x%in%tr$Species){
    return((tr$cont[which(tr$Species==x)]))    
  } else {
    return(as.numeric(NA))
  } } else { return(as.numeric(NA))}}

getTr2 <- function(x) {
  if(!is.na(x)){
    if(x %in% tr$Species){
      return(as.numeric(tr$cont[which(tr$Species==x)]))
    } else return(-1)
  } else return(as.numeric(NA))
}

ha <- readRDS("TraitsData/IUCN_PREDICTS_Habitats.rds")
ha2 <- readRDS("newIUCNHabitatinfo_pred.rds")
ha <- ha[,c(2:9,15)]
rs <- rowSums(ha,na.rm = T)
rs[rs==0] <- NA # If not occuring in at least one terrestrial habitat set to unknown
rs <- ifelse(rs==1,"habitat-specialist","habitat-generalist")
rs <- data.frame(Species=as.character(spec_lookup),ha=as.character(rs))
rs2 <- data.frame(Species=as.character(spec_lookup),ha=as.character(ifelse(ha2$ForestSpec=="Forest-specialist","Forest species","Non-forest species")))
getHa <- function(x) if(!is.na(x)) return(as.character(rs$ha[which(rs$Species==x)])) else return(NA)
getHa2 <- function(x) if(!is.na(x)) return(as.character(rs2$ha[which(rs2$Species==x)])) else return(NA)
brs <- readRDS("BirdLifeRangeSizes.rds")

res <- data.frame(SSS=unique(predictsdata$SSS),SpecINF_NarrowRanged=NA,SpecINF_AvgRange=NA,SpecINF_AvgRange2=NA,SpecINF_HabitatSpec=NA,SpecINF_ForestSpec=NA,SpecINF_IUCNStatus=NA,SpecINF_IUCNStatus2=NA)
for(study in unique(predictsdata$Source_ID)){
  print(paste("Processing",study))
  sub <- subset(predictsdata,Source_ID==study)
  sub$Best_guess_binomial <- droplevels(sub$Best_guess_binomial)
  sam <- dcast(sub,SSS~Best_guess_binomial,value.var = "Measurement",fun.aggregate = sum)
  rownames(sam) <- sam$SSS;sam$SSS <- NULL
  if(!length(names(sam))==1) { # No species inside?
    try( sam$Var.2 <- NULL ) # No idea whats inside? # Why does it create this column?
    sam <- repL(sam)
    # range
    comp <- apply(sam,1:2,getRange)
    val <- apply(comp,1,FUN=function(x) length(which(x=="narrow-ranged")))/ apply(comp,1,function(x) length(which(!is.na(x))))
    res$SpecINF_NarrowRanged[match(names(val),res$SSS)] <- val;rm(val)
    # AvgRange
    sp_c <- sam
    sp_c[sp_c==0] <- NA
    spc2 <- melt(as.matrix(sp_c),na.rm = T)
    spc2 <- merge(spc2,occ,by.x ="Var2",by.y="SpeciesID")
    d1 <- ddply(spc2,.(Var1),summarise,AvgRange=mean(log(NrCells),na.rm=T))
    res$SpecINF_AvgRange[match(d1$Var1,res$SSS)] <- d1$AvgRange
    rm(sp_c,spc2,d1)
    # AvgRange2
    sp_c <- sam
    sp_c[sp_c==0] <- NA
    spc2 <- melt(as.matrix(sp_c),na.rm = T)
    spc2 <- merge(spc2,brs,by.x ="Var2",by.y="SCINAME")
    d1 <- ddply(spc2,.(Var1),summarise,AvgRange=mean(log1p(DegOccupiedArea),na.rm=T))
    res$SpecINF_AvgRange2[match(d1$Var1,res$SSS)]  <- d1$AvgRange
    # Ha
    comp <- apply(sam,1:2,getHa)
    val <- apply(comp,1,FUN=function(x) length(which(x=="habitat-specialist")))/ apply(comp,1,function(x) length(which(!is.na(x))))
    res$SpecINF_HabitatSpec[match(names(val),res$SSS)] <- val;rm(val)        
    # HA-fa
    comp <- apply(sam,1:2,getHa2)
    val <- apply(comp,1,FUN=function(x) length(which(x=="Forest species")))/ apply(comp,1,function(x) length(which(!is.na(x))))
    res$SpecINF_ForestSpec[match(names(val),res$SSS)] <- val;rm(val)    
    # TR
    comp <- apply(sam,1:2,getTr)
    stopifnot(dim(sam)==dim(comp))
    try( val <- rowMeans(comp,na.rm = T) )
    if(exists("val")==F) val <- NA    
    res$SpecINF_IUCNStatus[match(names(val),res$SSS)] <- val;rm(val)
    # TR 2
    comp <- apply(sam,1:2,getTr2)
    stopifnot(dim(sam)==dim(comp))
    try( val <- rowMeans(comp,na.rm = T) )
    if(exists("val")==F) val <- NA else 
      val <- apply(comp,1,FUN=function(x) length(which(x >= 1))) / apply(comp,1,function(x) length(which(x >= 0)))
    res$SpecINF_IUCNStatus2[match(names(val),res$SSS)] <- val;rm(val)
    
    rm(sam,comp)
    }
}
# First load sites.pred
stopifnot(exists("sites.pred"))
stopifnot(length(which(!sites.pred$SSS%in%res$SSS)) ==0 )
# write.csv2(res,"CommunitywideTraitData.csv",row.names=F)
sites.pred <- plyr::join(sites.pred,res,by = "SSS") # ,all.x = T)
# Then kick out all nonsuitable studies and reclassify
rm(res)

#### Stats per Site ####
sites <- subset(sites,Grouping=="Aves")
library(marfunky)
# Subset PREDICTS only to Birds
sps <- subset(sites.pred,Grouping=="Birds")
field <- subset(sites,select=c("Transect","PREDICTS.LU","SpecINF_AvgRange","SpecINF_AvgRange2","SpecINF_HabitatSpec","SpecINF_ForestSpec","SpecINF_IUCNStatus2"),subset=Grouping=="Aves")
predicts <- subset(sps,select=c("PREDICTS.LU","SpecINF_AvgRange","SpecINF_AvgRange2","SpecINF_HabitatSpec","SpecINF_ForestSpec","SpecINF_IUCNStatus2"))
predicts <- predicts[which(!is.na(predicts$PREDICTS.LU)),]
predicts$PREDICTS.LU <- droplevels(predicts$PREDICTS.LU)
predicts$Transect <- "Africa-wide"
getTp <- function(lu,col){
  res = data.frame()
  for(l in lu){
    sub1 = subset(field,PREDICTS.LU==l)
    sub2 = subset(predicts,PREDICTS.LU==l)
    t = t.test(sub1[,col],sub2[,col])
    res <- rbind(res,data.frame(t=t$statistic,df=t$parameter,p=t$p.value))
  }
  return(res)
}
detachPackage("dplyr")
a = summarySEwithin(field,measurevar = "SpecINF_AvgRange2",withinvars = c("Transect","PREDICTS.LU"),na.rm = T);a$Group <- "Fieldwork"
a2 = summarySEwithin(predicts,measurevar = "SpecINF_AvgRange2",withinvars = c("Transect","PREDICTS.LU"),na.rm = T);a2$Group <- "Africa-wide"
a <- a[-grep("Urban",a$PREDICTS.LU),];a$PREDICTS.LU <- droplevels(a$PREDICTS.LU)# Exclude Urban
stopifnot(a$PREDICTS.LU==a2$PREDICTS.LU)
d <- rbind(a,a2);d <- cbind(d,getTp(d$PREDICTS.LU,"SpecINF_AvgRange2"));d$Type="African-wide range size\n (log-trans. area)"
names(d)[4:5] <- c("variable","variable2")

a= summarySEwithin(field,measurevar = "SpecINF_ForestSpec",withinvars = c("Transect","PREDICTS.LU"),na.rm = T);a$Group <- "Fieldwork"
a2= summarySEwithin(predicts,measurevar = "SpecINF_ForestSpec",withinvars = c("Transect","PREDICTS.LU"),na.rm = T);a2$Group <- "Broad-scale"
a <- a[-grep("Urban",a$PREDICTS.LU),];a$PREDICTS.LU <- droplevels(a$PREDICTS.LU)# Exclude Urban
stopifnot(a$PREDICTS.LU==a2$PREDICTS.LU)
d2 <- rbind(a,a2);d2 <- cbind(d2,getTp(d2$PREDICTS.LU,"SpecINF_ForestSpec"));d2$Type="Forest specialists\n (% of total)"
names(d2)[4:5] <- c("variable","variable2")

a = summarySEwithin(field,measurevar = "SpecINF_IUCNStatus2",withinvars = c("Transect","PREDICTS.LU"),na.rm = T);a$Group <- "Fieldwork"
a2 = summarySEwithin(predicts,measurevar = "SpecINF_IUCNStatus2",withinvars = c("Transect","PREDICTS.LU"),na.rm = T);a2$Group <- "Broad-scale"
a <- a[-grep("Urban",a$PREDICTS.LU),];a$PREDICTS.LU <- droplevels(a$PREDICTS.LU)# Exclude Urban
stopifnot(a$PREDICTS.LU==a2$PREDICTS.LU)
d3 <- rbind(a,a2);d3 <- cbind(d3,getTp(d3$PREDICTS.LU,"SpecINF_IUCNStatus2"));d3$Type="Threatened species\n (% of total)"
names(d3)[4:5] <- c("variable","variable2")

dd <- rbind(d,d2,d3)
dd$PREDICTS.LU <- ordered(dd$PREDICTS.LU,ord)
dd$Type <- ordered(dd$Type,c("African-wide range size\n (log-trans. area)","Forest specialists\n (% of total)","Threatened species\n (% of total)"))
dd <- dd[order(dd$PREDICTS.LU),]
#write.csv(dd,"/home/martin/Dokumente/Studium/Kopenhagen/Masters/000Thesis_writeup/Chapter1/TraitsTable.csv",row.names=F)

e <- as.data.frame(dd)
e$pval <- ifelse(e$p<0.001,"*","")
ee <- melt(e,id.vars = c("Transect","PREDICTS.LU","Group","Type","se","pval"),measure.vars = c("variable") )
ee$Type <- as.factor(ee$Type)

### Figure -
rect_lab <- c("PV","SV","PL","CL")#,"UR")
rect_labf <- c("Primary\n vegetation","Secondary\n vegetation","Plantation\n forest","Cropland")
rect_colours<-c("#00AE00","#94BD5E", "#006B6B", "#E6E64C")#,"#808080")
xmin = c(0.5,1.5,2.5,3.5);xmax=c(1.5,2.5,3.5,4.5)

g <- ggplot(ee,aes(x=PREDICTS.LU,y=value,fill=Transect)) + theme_few(base_size = 16)
g <- g + geom_bar(stat="identity",position="dodge",width=.75,alpha=.6)
g <- g + geom_errorbar(aes(ymin=value-se, ymax=value+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))
g <- g + facet_wrap(~Type,ncol = 1,scales = "free_y") + theme(strip.text.x = element_text(size=26))#,strip.background = element_rect(colour="black"))
#g <- g + scale_y_continuous(expand=c(0,0),limits=c(0,.55),breaks=pretty_breaks())
#g <- g + geom_text(aes(y=.5,label=pval),size=10)
g <- g + scale_x_discrete(labels=rect_labf) + theme(axis.title.x=element_text(size=26))
g <- g + theme(panel.margin = unit(0.5, "lines"))
g <- g + xlab("") + ylab("Community wide average") # No labels  
g <- g + scale_fill_manual(values=c("red","blue","black"),guide = guide_legend(reverse = T,direction = "horizontal", title.position = "top",title="",ncol=1,
                                              label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                              label.theme = element_text(angle = 90,colour="black",size=24))) 
g <- g + theme(axis.text.x = element_text(size = rel(1.2), colour = 'black')) +theme(axis.text.y = element_text(size = 24, colour = 'black'))
g <- g + theme(axis.ticks.x=element_blank()) # Remove legend and x-axis stuff
g <- g + theme(legend.position = "none")
g
ggsave("Figure5.png",plot= g,scale = 1.1,units = "mm",dpi=400)
ggsave("../000Thesis_writeup/WriteUpLatex/gfx/Figure5.pdf",plot=g,scale=1.1,dpi=400)

# Transparent
t <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "transparent",colour = NA),
           axis.text = element_text(colour="black"),axis.text.y =element_text(colour="black"), axis.line = element_line(colour="black"),axis.title = element_text(colour="black"),
           legend.key = element_rect(fill=NA,colour = NA),legend.background = element_rect(fill = NA,colour=NA),legend.text = element_text(color="black"),
           title = element_text(colour="black"),
           strip.background = element_rect(fill=NA,color=NA), strip.text = element_text(color="black"),
           panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
           text = element_text(color="black",size=24)
) 
ggsave("Figure5_trans.png",plot=g+t,scale = 1.1,units = "mm",dpi=400,bg="transparent")


#### Missing trait data ####
# Field
occ <- readRDS("TraitsData/GBIF_OCC-Fielddata.rds") #Range
occ$NrCells[which(occ$NrCells==0)] <- NA
a1 <- 1-(length(which(is.na(occ$NrCells))) / length(occ$NrCells))

f1 <- read.csv("TraitsData/IUCN_field_data.csv")
a2 <- 1-(length(which(f1$Status=="NE")) / length(f1$Status))

f <- readRDS("newIUCNHabitatinfo.rds")
1-(length(which(is.na(f$ForestSpec))) / length(f$ForestSpec))


## PREDICTS
occ <- readRDS("TraitsData/GBIF_OCC-PREDICTS.rds")
occ$NrCells[which(occ$NrCells==0)] <- NA
b1 <- 1-(length(which(is.na(occ$NrCells))) / length(occ$NrCells))

f <- readRDS("newIUCNHabitatinfo_pred.rds")
b2 <- 1-(length(which(is.na(f$ForestSpec))) / length(f$ForestSpec))

f1 <- read.csv("TraitsData/IUCN_PREDICTS_Threats.csv",header=T)
1-(length(which(f1$Status=="NE")) / length(f1$Status))

d <- data.frame(matrix(list(a1,a2,b1,b2),nrow = 2))
colnames(d) <- c("Field","PREDICTS")
rownames(d) <- c("GBIF","IUCN")


theme_set(theme_classic(base_size=16,base_family = "sans"))
g <- ggplot(ee,aes(x=PREDICTS.LU,y=value,fill=Transect))
g <- g + geom_bar(stat="identity",position="dodge",width=.75,alpha=.6)
g <- g + geom_errorbar(aes(ymin=value-se, ymax=value+se),
                       width=.2,                    # Width of the error bars
                       position=position_dodge(.9))
g <- g + facet_wrap(~Type,ncol = 1,scales = "free_y") + theme(strip.text.x = element_text(size=14))#,strip.background = element_rect(colour="black"))
#g <- g + scale_y_continuous(expand=c(0,0),limits=c(0,.55),breaks=pretty_breaks())
#g <- g + geom_text(aes(y=.5,label=pval),size=10)
#g <- g + scale_x_discrete(labels=rect_lab)
g <- g + scale_x_discrete(labels=rect_lab) + theme(axis.title.x=element_text(size=20))
g <- g + theme(panel.margin = unit(0.5, "lines"))
g <- g + xlab("") + ylab("Community wide average") # No labels  
g <- g + scale_fill_brewer(palette = "Dark2",guide = guide_legend(reverse = T,direction = "horizontal", title.position = "top",title="",ncol=1,
                                                                  label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                                                  label.theme = element_text(angle = 90,colour="white"))) 
g <- g + theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) # Remove legend and x-axis stuff
g
ggsave("OutputComparison/AverageCommunityTraits.png",plot= g + scale_fill_grey(),scale = 1.1,units = "mm",dpi=400)
ggsave("../000Thesis_writeup/WriteUpLatex/gfx/Figure5.pdf",plot=g,scale=1.1,dpi=400)
