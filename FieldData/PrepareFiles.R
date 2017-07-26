library(gdata)
library(BiodiversityR)
library(vegan)
dir = getwd()
#### Taita ####
print("...Reading in my own Data...")
setwd("FieldData/")
taita_1 <- read.xls("Species.xls",sheet=2,perl="C:/Perl64/bin/perl.exe")
taita_2 <- read.xls("Species.xls",sheet=3,perl="C:/Perl64/bin/perl.exe")
taita_1[is.na(taita_1)] <- 0 # Replace NA with 0
taita_2[is.na(taita_2)] <- 0 # Replace NA with 0
rownames(taita_1) <- taita_1[,1]
rownames(taita_2) <- taita_2[,1]
taita_1 <- taita_1[,-1] # Kick out species column
taita_2 <- taita_2[,-1] # Kick out species column

taita <- taita_1 + taita_2 # Both sampling dates together
#### Kili ####
kili_1 <- read.xls("Species.xls",sheet=4,perl="C:/Perl64/bin/perl.exe")
kili_2 <- read.xls("Species.xls",sheet=5,perl="C:/Perl64/bin/perl.exe")
kili_1[is.na(kili_1)] <- 0 # Replace NA with 0
kili_2[is.na(kili_2)] <- 0 # Replace NA with 0
rownames(kili_1) <- kili_1[,1]
rownames(kili_2) <- kili_2[,1]
kili_1 <- kili_1[,-1] # Kick out species column
kili_2 <- kili_2[,-1] # Kick out species column

kili <- kili_1 + kili_2 # Both sampling dates together

# Merge both
rm(kili_1,kili_2,taita_1,taita_2)
taita$Spec <- row.names(taita)
kili$Spec <- row.names(kili)
species <- merge(taita,kili,by="Spec",all=T)
species[is.na(species)] <- 0 # Replace NA with 0
row.names(species) <- species$Spec # set final row.names
# kick out names
species$Spec <- NULL
taita$Spec <- NULL
kili$Spec <- NULL
species <- t(species)
rm(kili,taita)



# Site information
print("...Generate site-level information...")
sites <- read.xls("SiteData.xls",sheet=1,perl="C:/Perl64/bin/perl.exe")

sites$Landuse <- as.factor(sites$Landuse)
sites$PREDICTS.LU <-as.factor(sites$PREDICTS.LU)
sites$PREDICTS.LUI <- as.factor(sites$PREDICTS.LUI)
sites$abund[match(rownames(species),sites$SiteName)] <- apply(species,1,function(x) sum(x,na.rm = T))
sites$spec[match(rownames(species),sites$SiteName)] <- rowSums(species > 0,na.rm = T)
sites$shan <-  diversitycomp(species,sites,factor1="SiteName",index="Shannon",method="all")[,2][match(sites$SiteName,names(diversitycomp(species,sites,factor1="SiteName",index="Shannon",method="all")[,2]))]
sites$simp <-  diversitycomp(species,sites,factor1="SiteName",index="inverseSimpson",method="all")[,2][match(sites$SiteName,names(diversitycomp(species,sites,factor1="SiteName",index="inverseSimpson",method="all")[,2]))]

sites$logabund <- log1p(sites$abund)
sites$logsimp <- log1p(sites$simp)


# Get correct start/end date
sites$Sample_start_earliest <- paste0(formatC(sites$DateFirst,width = 5,flag="0"),".2014")
sites$Sample_end_latest <- paste0(formatC(sites$DateSecond,width = 4,flag="0"),".2014")

sites$Sample_start_earliest <-  as.Date(sites$Sample_start_earliest,"%d.%m.%Y")
sites$Sample_end_latest <-  as.Date(sites$Sample_end_latest,"%d.%m.%Y")

sites$Grouping <- "Aves"

#### Extract data for Dickens ####
#dick <- sites[which(sites$Authority=="DickensO"),] # only dickens sites 
#plot(dick$abund~dick$Elev,col=dick$Transect)
#summary(lm((dick$abund)~dick$Elev))
# rm(dick)
# 
#dick <- sites[which(sites$Authority=="DickensO"),] # only dickens sites 
#dick <- species[which(rownames(species)%in%dick$SiteName[which(dick$Authority=="DickensO")]),]
#dick <- dick[order(rownames(dick)),]
#dick <- dick[,colSums(dick)!=0]# kick empty species out
#dick <- ifelse(dick == 0,0,1) # make Occurence grid
#sp <- read.xls("Species.xls",1)
#colnames(dick) <- sp$Scientific.Name[match(colnames(dick),sp$Common.Name)]
#write.csv2(dick,"MJ_DickensSites.csv")

# Load in Dickens data
print("...Reading in Dickens Data...")
plantdata <- read.xls("Dickens_Taita_Kilimanjaro_Stock_Density.xls",1,perl="C:/Perl64/bin/perl.exe")
plantdata <- plantdata[,-ncol(plantdata)];plantdata$PlotNumb <- as.character(plantdata$PlotNumb)
rownames(plantdata) <- plantdata$PlotNumb;plantdata$PlotNumb <- NULL

print("...Generate site-level information for dickens data...")
sites.dickens <- read.xls("SiteData.xls",sheet=1,perl="C:/Perl64/bin/perl.exe");sites.dickens <- subset(sites.dickens,Authority=="DickensO")
sites.dickens$ID <- seq(nrow(sites)+1,nrow(sites)+nrow(sites.dickens))
sites.dickens$Landuse <- as.factor(sites.dickens$Landuse)
sites.dickens$PREDICTS.LU <-as.factor(sites.dickens$PREDICTS.LU)
sites.dickens$PREDICTS.LUI <- as.factor(sites.dickens$PREDICTS.LUI)

sites.dickens$abund[match(rownames(plantdata),sites.dickens$SiteName)] <- apply(plantdata,1,function(x) sum(x,na.rm = T))
sites.dickens$spec[match(rownames(plantdata),sites.dickens$SiteName)] <- rowSums(plantdata > 0,na.rm = T)
sites.dickens$shan <-  diversitycomp(plantdata,sites.dickens,factor1="SiteName",index="Shannon",method="all")[,2][match(sites.dickens$SiteName,names(diversitycomp(plantdata,sites.dickens,factor1="SiteName",index="Shannon",method="all")[,2]))]
sites.dickens$simp <-  diversitycomp(plantdata,sites.dickens,factor1="SiteName",index="inverseSimpson",method="all")[,2][match(sites.dickens$SiteName,names(diversitycomp(plantdata,sites.dickens,factor1="SiteName",index="inverseSimpson",method="all")[,2]))]
sites.dickens$logabund <- log1p(sites.dickens$abund)
sites.dickens$logsimp <- log1p(sites.dickens$simp)

sites.dickens$Sample_start_earliest <- NA
sites.dickens$Sample_end_latest <- NA
sites.dickens$Sample_start_earliest[sites.dickens$Transect=="Taita"] <- as.Date(strptime("01-09-2011","%d-%m-%Y"),origin = "1900-01-01")
sites.dickens$Sample_start_earliest[sites.dickens$Transect=="Kilimanjaro"] <- as.Date(strptime("01-04-2012","%d-%m-%Y"),origin = "1900-01-01")
sites.dickens$Sample_end_latest[sites.dickens$Transect=="Taita"] <- as.Date(strptime("01-10-2011","%d-%m-%Y"),origin = "1900-01-01")
sites.dickens$Sample_end_latest[sites.dickens$Transect=="Kilimanjaro"] <- as.Date(strptime("01-05-2012","%d-%m-%Y"),origin = "1900-01-01")

#I believe that internally, dates are stored as numeric seconds since the start of Unix time (1st January 1970)
sites.dickens$Sample_start_earliest <- as.Date(sites.dickens$Sample_start_earliest, origin = "1970-01-01")
sites.dickens$Sample_end_latest<- as.Date(sites.dickens$Sample_end_latest, origin = "1970-01-01")

#sites.dickens$Sample_start_earliest <- paste0(formatC(sites.dickens$DateFirst,width = 5,flag="0"),".2014")
#sites.dickens$Sample_end_latest <- paste0(formatC(sites.dickens$DateSecond,width = 4,flag="0"),".2014")
#sites.dickens$Sample_start_earliest <-  as.Date(sites.dickens$Sample_start_earliest,"%d.%m.%Y")
#sites.dickens$Sample_end_latest <-  as.Date(sites.dickens$Sample_end_latest,"%d.%m.%Y")

sites.dickens$Grouping <- "Plantae"
# FINAL
print("Merge both and cleanup")
sites <- rbind(sites,sites.dickens)
rm(sites.dickens)
species.plants <- as.matrix(plantdata)

setwd(dir)
rm(dir)