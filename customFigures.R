#### Preperation and standard Package Loading ####
#rm(list=ls())
par.ori <- par(no.readonly=T)
library(ggplot2) # advanced fancy plotting
library(scales)
#library(extrafont)
library(rgdal)
library(raster)
library(rasterVis)
library(rgeos)
library(gridExtra)
#library(ggmap)

sp <- sites.pred
coordinates(sp) <- ~Longitude + Latitude
proj4string(sp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
writeOGR(sp,"GIS","PREDICTS_Africa",driver = "ESRI Shapefile")

dir.create("presentation_figures")
# To be run for individual Figures outside analysis stuff
transect_t <- readOGR("../GIS_Field/","Taita_transectbuffer_1km_Arc1960_37s")
transect_k <- readOGR("../GIS_Field/","Kilimanjaro_transect_buffer1km_arc1960")
# Transform to WGS84
ll <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
transect_t <- spTransform(transect_t,ll)
transect_k <- spTransform(transect_k,ll)

kilitaita <- readOGR("../GIS_Field/","Extent_nice_big") #Just for extent rectangle
ex <- extent(kilitaita)
cropNice <- readOGR("../GIS_Field/","Extent_nice")

source("../FieldData/PrepareFiles.R")
source("GIS_Loading.R")
sites <- field_data_load(sites)
sites$PREDICTS.LU <- ordered(sites$PREDICTS.LU,c("Primary Vegetation", "Secondary Vegetation","Plantation forest","Pasture", "Cropland","Urban"))

# Get Background ForestCover from Hansen
fc <- raster("../GIS/HANSEN/TC2000_KiliTaitaClip.tif")
fc_t <- crop(fc,subset(cropNice,loc=="Taita"))
fc_k <- crop(fc,subset(cropNice,loc=="Kili"))

dem_t <- raster("../GIS/DEM/TaitaDEM_Pub.tif")
hill_t <- raster("../GIS/DEM/Taita_CropNice_hillshade.tif")
dem_t <- crop(dem_t,subset(kilitaita,loc=="Taita"));hill_t <- crop(hill_t,subset(kilitaita,loc=="Taita"));
dem_k <- raster("../GIS/DEM/KiliDEM_Pub.tif")
hill_k <- raster("../GIS/DEM/Kili_CropNice_hillshade.tif")
dem_k <- crop(dem_k,subset(kilitaita,loc=="Kili"));hill_k <- crop(hill_k,subset(kilitaita,loc=="Kili"));


#aa <- SpatialPoints(cbind(sites$Longitude,sites$Latitude),proj4string = ll)
#dem_t <- crop(dem_t,bbox(aa));hill_t <- crop(hill_t,bbox(aa))
#dem_k <- crop(dem_k,bbox(aa));hill_k <- crop(hill_k,bbox(aa))

#slope <- terrain(dem_t, opt='slope')
#aspect <- terrain(dem_t, opt='aspect')
#hill_t <- hillShade(slope, aspect, 40, 270)

#slope <- terrain(dem_k, opt='slope')
#aspect <- terrain(dem_k, opt='aspect')
#hill_k <- hillShade(slope, aspect, 40, 270)

# Map of study area
# Show both transects near other + Forest cover
# Transect maps seperated by line
# Inspired here
# http://r-nold.blogspot.dk/2014/06/creating-inset-map-with-ggplot2.html
KEN <-raster::getData("GADM", country="KEN", level=0) 
TZ <-raster::getData("GADM", country="TZA", level=0) 
row.names(TZ) <- "2"
EA <- sp::rbind.SpatialPolygonsDataFrame(KEN,TZ,makeUniqueIDs=T) # Combine them

# Extent rectangle for inset map
pol<-data.frame(xmin=ex@xmin,xmax=ex@xmax ,ymin=ex@ymin ,ymax=ex@ymax)

#### Figure 1  - Map building ####
#ggthemr("earth",layout = "clean",type = "outer",spacing = 0)
theme_set(theme_classic(base_size=12,base_family = "sans"))

rcols <- c("#00AE00","#94BD5E", "#006B6B", "#E6E64C","#FF8000")
# Main Map - Kili
gm <- ggmap( get_googlemap(zoom=12,maptype="roadmap",center=cbind(37.4006,-3.3651),
                           scale=2))

tk = fortify(transect_k)
hk <- as.data.frame(rasterToPoints(hill_k));names(hk) <- c("Longitude", "Latitude", "value")
dk <- as.data.frame(rasterToPoints(dem_k));names(dk) <- c("Longitude", "Latitude", "value")
#pk <- gm
pk <- ggplot()
# Hillshade + Altitude
pk <- pk + geom_raster(data=dk,aes(y=Latitude,x=Longitude,fill=value),alpha=0.75)
pk <- pk + geom_raster(data=hk,aes(y=Latitude,x=Longitude,alpha=1-value),fill="gray20") + scale_alpha(guide=FALSE,range = c(0,1.00))
pk <- pk + scale_fill_gradientn(name="Altitude",colours = grey.colors(100),guide="none")
pk <- pk + coord_equal()
pk <- pk + geom_polygon(data=tk, aes(long ,lat, group=group), colour="grey10",fill="grey80",alpha=.5)
pk <- pk + geom_point(data=sites[which(sites$Transect=="Kilimanjaro"),],
                      aes(x=Long,y=Lat,color=PREDICTS.LU),size=4,shape=15,alpha=.9)
pk <- pk + xlab("Longitude (°)") + ylab("")
pk <- pk + scale_y_continuous(breaks=pretty_breaks())+scale_x_continuous(breaks=pretty_breaks())
#pk <- pk + ggtitle("Kirua Vunjo Transect (Kilimanjaro)")
pk <- pk + scale_color_manual(values = rcols,guide=guide_legend(title="Land use",hjust=.5))
pk <- pk + theme(legend.margin=unit(0,"lines"), legend.box="vertical",
                 legend.key.size=unit(1,"lines"), legend.text.align=0,
                 legend.title.align=0)
pk <- pk + theme(axis.text=element_text(size=16), axis.title = element_text(size=18)) # adjust text size
pk
#pk <- pk + theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),panel.background=element_blank(),panel.grid=element_blank())
#ggsave(filename="presentation_figures/KiliMap.png",plot=pk+theme(legend.position="none"),dpi=400)

# Main Map - Taita
gt <- ggmap( get_googlemap(zoom=12,maptype="terrain",center=coordinates(gCentroid(SpatialPoints(cbind(sites$Long[which(sites$Transect=="Taita")],sites$Lat[which(sites$Transect=="Taita")])))),
                           scale=2))

tt = fortify(transect_t)
ht <- as.data.frame(rasterToPoints(hill_t));names(ht) <- c("Longitude", "Latitude", "value")
dt <- as.data.frame(rasterToPoints(dem_t));names(dt) <- c("Longitude", "Latitude", "value")

#pt <- gt
pt <- ggplot()
pt <- pt + geom_raster(data=dt,aes(y=Latitude,x=Longitude,fill=value),alpha=0.75)
pt <- pt + geom_raster(data=ht,aes(y=Latitude,x=Longitude,alpha=1-value),fill="gray20") + scale_alpha(guide=FALSE,range = c(0,1.00))
pt <- pt + scale_fill_gradientn(name="Altitude",colours = grey.colors(100),guide="none")
pt <- pt +  coord_equal()
pt <- pt + geom_polygon(data=tt, aes(long ,lat, group=group), colour="grey10",fill="grey80",alpha=.5)
pt <- pt + geom_point(data=sites[which(sites$Transect=="Taita"),],
                      aes(x=Long,y=Lat,color=PREDICTS.LU),size=4,shape=15,alpha=.9)

pt <- pt+ xlab("Longitude (°)") + ylab("Latitude (°)")
pt <- pt + scale_x_continuous(breaks=pretty_breaks(4))
#pt <- pt + ggtitle("Wundanyi Transect (Taita Hills)")
pt <- pt + scale_color_manual(values = rcols,guide=guide_legend(title="Land use",hjust=.5))
pt <- pt + theme(axis.text=element_text(size=16),axis.title = element_text(size=18)) # adjust text size
pt
#ggsave(filename="presentation_figures/TaitaMap.png",plot=pt+theme(legend.position="none"),dpi=400)


# Inset Map
pi <- ggplot()
pi <- pi + geom_polygon(data=EA, aes(long,lat,group=group),colour="black",size=1.2,fill="#fff7bc")
pi <- pi + coord_equal()+labs(x=NULL,y=NULL)
pi <- pi + geom_rect(data = pol, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha=0.5,fill="red", colour="red", size = 1.2, linetype=1)
pi <- pi +theme(axis.text.x =element_blank(),axis.text.y= element_blank(), axis.ticks=element_blank(),axis.title.x =element_blank(),
                axis.title.y= element_blank(),axis.line=element_blank(),axis.line=element_blank())
#ggsave(filename="presentation_figures/InsetMap.png",plot=pi,width = 12,height=12,unit="in",dpi=400)

# Transparent
t <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "transparent",colour = NA),
           axis.text = element_text(colour="white"),axis.text.x =element_text(colour="white"),axis.text.y =element_text(colour="white"),axis.ticks = element_line(colour = "white"), axis.line = element_line(colour="white"),axis.title = element_text(colour="white"),
           legend.key = element_rect(fill=NA,colour = NA),legend.background = element_rect(fill = NA,colour=NA),legend.text = element_text(colour="white"),
           panel.grid.minor = element_blank(), title = element_text(colour="white")
)


ggsave(filename="../001Thesis_Defense/assets/img/InsetMap_trans.png",plot=pi+t,scale=1.1,dpi=400,bg = "transparent")
ggsave(filename="../001Thesis_Defense/assets/img/Kilimap_trans.png",plot=pk+t,scale=1.1,dpi=400,bg = "transparent")
ggsave(filename="../001Thesis_Defense/assets/img/Taitamap_trans.png",plot=pt+t,scale=1.1,dpi=400,bg = "transparent")

# Export gridextra
# On one page with combined Insert + legend
# library(marfunky)
# mylegend<-g_legend(pk)
# png("presentation_figures/legend.png",width=16,height=16,units="in",res=400)
# grid.arrange(mylegend)
# dev.off()
ggsave("../Conferences/Conservation Symposium/taita.png",pt + theme(legend.position="none",plot.margin=unit(c(1,1,0,1), "cm")))
ggsave("../Conferences/Conservation Symposium/kili.png",pk + theme(legend.position="none",plot.margin=unit(c(-0,1,1,1), "cm")))

# Now combine the shit with grobs
svg("presentation_figures/CombinedMaps.svg",width=8,height=6)
#png("../Conferences/Conservation Symposium/FieldWork.png",width = 12,height=12,res = 400,units = "in",bg="black")
#svg("presentation_figures/CombinedMaps.svg",width = 18,height=18)
grid.arrange(clip = T,arrangeGrob(
  pt + theme(legend.position="none",plot.margin=unit(c(1,1,0,1), "cm")),
  pk + theme(legend.position="none",plot.margin=unit(c(-0,1,1,1), "cm")),
  nrow=2,heights=c(1,1))
#  ,
#             arrangeGrob(clip = T, mylegend)
                        , ncol=2,nrow=1,heights=c(1,2,3,3)
)
dev.off()

# Do the rest in inkscape, shezzzz

# Export
png(file="presentation_figures/overviewMap.png",w=3000,h=1500, res=400)
grid.newpage()

v1<-viewport(width = 0.5, height = 0.5, x = 0, y = 0.5) #plot area for the main map
v2<-viewport(width = 0.5, height = 0.5, x = 0.5, y = 0.5) #plot area for the main map
vi<-viewport(width = 0.3, height = 0.3, x = 0.86, y = 0.28) #plot area for the inset map

print(pk,vp=v1)
print(pt,vp=v2)
print(pi,vp=vi)

dev.off()

#### PREDICTS STUDIES SI TABLE 1####
library(stargazer)
library(stringr)
tab <- read.csv("../000Thesis_writeup/PaperDraft/Appendix/SI_Table1.csv")
tab$DOI <- as.character(tab$DOI)
tab$DOI <- ifelse(str_length(tab$DOI)>0,paste0("doi{",tab$DOI,"}"),"")
tab$Year_of_publication <- as.factor(tab$Year_of_publication)
tab$Journal_title
tab$Journal_title <- str_replace_all(tab$Journal_title,"&","\\&")
names(tab) <- c("First author surname","Publication year","Publication title","Journal","DOI","Sites (N)","Land use classes (N)","Investigated Taxon","Species richness (N)")
tab$DOI <- NULL
stargazer(tab,type = "latex",summary = F,header=F,font.size = "tiny",rownames = F,label = "tab_app:Table1",column.sep.width = "2pt",
          float.env="sidewaystable")

#### Predicts data general map ####
library(maps)
library(rgeos)
africa <- readOGR("GIS","Country_Africa")
africa.map = ggplot2::fortify(africa, region="COUNTRY_NA")
# Map including flipped stacked barplot of groupings
sites.pred$Lat.r <- round(sites.pred$Latitude)
xcol <- brewer.pal(7,"Dark2")

library(ggthemr)
ggthemr("pale",layout = "minimal",type = "outer",spacing = 0)
#theme_set(theme_classic(base_size=12,base_family = "sans"))
hr <- ggplot(sites.pred,aes(x=Lat.r,group=Grouping,fill=Grouping))
hr <- hr + geom_bar(position="stack")+coord_flip() + theme(legend.position="none",axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                                           axis.title.x=element_blank(), axis.title.y=element_blank())
hr <- hr + scale_x_continuous(breaks=seq(-30,40,10))
hr <- hr + scale_fill_manual(values = xcol)
hr <- hr + labs(x="",y="") 
hr <- hr + theme(plot.margin=unit(c(1,1,1,-.5), "cm")) # Remove margins

library(dplyr)
kt <- africa.map %>% filter(id %in% c("Kenya","Tanzania, United Rep of"))
g <- ggplot(africa.map, aes(x = long, y = lat, group = group))
g <- g + geom_polygon(colour = "grey80", size = 0.5, fill = "white", aes(group = group))
#g <- g + geom_polygon(data=kt,colour = "black", size = 0.8, fill = "#fff7bc",alpha=0.75, aes(group = group))
g <- g + geom_point(data = sites.pred, aes(x = Longitude, y = Latitude, fill = Grouping, group = NULL), shape = 21, colour = "black", size = 3)
g <- g + scale_fill_manual(values = xcol)
g <- g + xlab("Longitude (°)") + ylab("Latitude (°)")
g <- g + scale_y_continuous(expand=c(0,1),limits=c(-38,40),breaks=seq(-30,40,10))
# annotate a rectangle
#a <- as.vector(bbox(SpatialPoints(cbind(sites$Long,sites$Lat))))
#g <- g + geom_rect(xmin=37,ymin=-5,xmax=39,ymax=0,color="red",fill=NA)
g <- g + theme(legend.position=c(.2,.2),legend.title=element_blank(),legend.background=element_blank(),legend.key.height=unit(c(.5),"cm"))
g <- g + theme(plot.margin=unit(c(1,0,1,1), "cm"), # Remove margins
               panel.border = element_rect(colour = "black",fill = NA,size=1.5),
               panel.grid.minor=element_line(size=.1,colour="grey")) 

gg <- arrangeGrob(g, hr, ncol=2, nrow=1, widths=c(4,1), heights=c(2, 1))

ggsave(filename="OutputComparison/PREDICTS_Africa.png",plot=gg,scale=1.1,dpi=400)
# Final thesis
ggsave(filename="../000Thesis_writeup/WriteUpLatex/gfx/SIFigure1.pdf",plot=gg,dpi=400)

# Plot transparent
t <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "transparent",colour = NA),
           axis.text = element_text(colour="white"),axis.text.x =element_text(colour="white"),axis.text.y =element_text(colour="white"),axis.ticks = element_line(colour = "white"), axis.line = element_line(colour="white"),axis.title = element_text(colour="white"),
           legend.key = element_rect(fill="white",colour = "white"),legend.background = element_rect(fill = "white",colour="white"),
           panel.grid.minor = element_blank()
           )
gg <- arrangeGrob(g + t,
                  hr + t + theme(axis.text=element_blank()),
                  ncol=2, nrow=1, widths=c(4,1), heights=c(2, 1))

ggsave(filename="../001Thesis_Defense/assets/img/PREDICTS_Africa_trans.png",plot=gg,scale=1.1,dpi=400,bg = "transparent")



#### SpaMM Richness and Abundance distribution ####
library(rgeos);library(spaMM)
# spaMM package for correlated error terms
names(sub.sites)
sub.sites$logpop <- log(sub.sites$pop) # Logtransform Population
sub.s <- subset(sub.sites,Transect=="Kilimanjaro")
coordinates(sub.s) <- ~ Lat+Long; proj4string(sub.s) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
# Only Points on Transect
sub.s <- sub.s[which(gWithin(sub.s,transect_k,byid = T)),]
sub.s <- as.data.frame(sub.s);names(sub.s)[1:2] <- c("Lat","Long")

ga_full <- glm(spec~PREDICTS.LU+Elev+I(Elev^2)+logpop+I(logpop^2)+meanEVI+I(meanEVI^2),data=sub.s,family=poisson)
up <- stepAIC(ga_full,direction="back",trace=T)

fst <- as.formula(update.formula(up$formula,"~. + Matern(1|Lat+Long)"))
sfit_lu <- corrHLfit(fst,data=sub.s,family=gaussian(),HLmethod="ML")


png("OutputDiversityAnalyis/SPAMM_Modelling_Kili_Richness.png",width = 20,height=9,units = "in",res = 400)
filled.mapMM(sfit_lu,nlevels=30,plot.axes={axis(1);axis(2)},color.palette = function(n){spaMM.colors(n,redshift=1/2)},
             plot.title=title(main="Sampling Kilimanjaro",
                              xlab="longitude",ylab="latitude"))
dev.off()

#### NOTE: Insert bird study points + NR studies on OVERLAP #####

#### PREDICTS scatterplot ####
source("study_specific_models.R")
sem <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
#keep urban
sub.sites <- sites#subset(sites,PREDICTS.LU!="Urban")
sub.sites$LogAbund <- sub.sites$logabund
sub.sites$LogSimpson <- sub.sites$logsimp
sub.sites$Species_richness <- sub.sites$spec

#sub.sites.pred <- subset(sites.pred,PREDICTS.LU!=c("Pasture")) # Pasture out
sub.sites.pred <- sites.pred
sub.sites.pred$PREDICTS.LU <-droplevels(sub.sites.pred$PREDICTS.LU)
sub.sites.pred$PREDICTS.LUI[which(sub.sites.pred$PREDICTS.LUI=="Cannot decide")] <- NA
sub.sites.pred$PREDICTS.LUI <-droplevels(sub.sites.pred$PREDICTS.LUI)

sub.sites.pred$LogAbund <-sub.sites.pred$logabund
sub.sites.pred$ChaoRint <-sub.sites.pred$ChaoR
sub.sites.pred$LandUse <- sub.sites.pred$PREDICTS.LU
sub.sites.pred$UseIntensity <- sub.sites.pred$PREDICTS.LUI
sub.sites.pred$UI <- as.factor(paste(sub.sites.pred$LandUse,sub.sites.pred$UseIntensity))

out <- study_specific_models(sub.sites.pred,glmer.sp.richness=T,outDir="/tmp",graphics.format="png",graphics.res=400)
# Mark bird studies in consitency
out$Grouping <- NA
birds <- sub.sites.pred$SS[which(sub.sites.pred$Grouping=="Birds")]
#plants <- sub.sites.pred$SS[which(sub.sites.pred$Grouping=="Plants")]
out$Grouping <- ifelse(out$studyID%in%birds,"Birds",NA)
#out$Grouping <- ifelse(out$studyID%in%plants,"Plants",NA)
rm(birds)

# FiT GLM for Fielddata
formatData <- function(afit,resp_var="LogAbund",gr="Taita"){
  require(arm)
  land_use <- ord[c(-4)]
  land_use <- ordered(land_use,levels=(ord[c(-4)]))
  estimate <- coef(afit);estimate[1] = 0
  se <- se.coef(afit)
  trans_est<-(exp(estimate))*100
  upper_est<-(exp(estimate+se))*100
  lower_est<-(exp(estimate-se))*100
  lu <- sort(unique(afit$data$PREDICTS.LU))
  d1 <- data.frame(lu,estimate,se,lower_est,trans_est,upper_est)
  d1$studyID <- gr
  d1$resp_var <- resp_var
  return(d1)
}
# tAITA
afit <- glm(LogAbund~PREDICTS.LU,data=subset(sub.sites,Transect=="Taita"&Grouping=="Aves"),family="gaussian")
sfit <- glm(Species_richness~PREDICTS.LU,data=subset(sub.sites,Transect=="Taita"&Grouping=="Aves"),family=poisson(link=log))
a1 <- formatData(afit,"LogAbund","Taita");a1$Grouping <- "Birds"
b1 <- formatData(sfit,"Species_richness","Taita");b1$Grouping <- "Birds"
#afit <- glm(LogAbund~PREDICTS.LU,data=subset(sub.sites,Transect=="Taita"&Grouping=="Plantae"),family="gaussian")
#sfit <- glm(Species_richness~PREDICTS.LU,data=subset(sub.sites,Transect=="Taita"&Grouping=="Plantae"),family=poisson(link=log))
#a2 <- formatData(afit,"LogAbund","Taita");a2$Grouping <- "Plants"
#b2 <- formatData(sfit,"Species_richness","Taita");b2$Grouping <- "Plants"
# Kili
afit <- glm(LogAbund~PREDICTS.LU,data=subset(sub.sites,Transect=="Kilimanjaro"&Grouping=="Aves"),family="gaussian")
sfit <- glm(Species_richness~PREDICTS.LU,data=subset(sub.sites,Transect=="Kilimanjaro"&Grouping=="Aves"),family=poisson(link=log))
c1 <- formatData(afit,"LogAbund","Kilimanjaro");c1$Grouping <- "Birds"
d1 <- formatData(sfit,"Species_richness","Kilimanjaro");d1$Grouping <- "Birds"
#afit <- glm(LogAbund~PREDICTS.LU,data=subset(sub.sites,Transect=="Kilimanjaro"&Grouping=="Plantae"),family="gaussian")
#sfit <- glm(Species_richness~PREDICTS.LU,data=subset(sub.sites,Transect=="Kilimanjaro"&Grouping=="Plantae"),family=poisson(link=log))
#c2 <- formatData(afit,"LogAbund","Kilimanjaro");c2$Grouping <- "Plants"
#d2 <- formatData(sfit,"Species_richness","Kilimanjaro");d2$Grouping <- "Plants"

d1b <- rbind(a1,b1);d1b$Type<-"Fieldwork";rm(a1,b1,a2,b2)
d1p <- rbind(c1,d1);d1p$Type<-"Fieldwork";rm(c1,d1,c2,d2)
d1 <- rbind(d1b,d1p);d1$land_use <- d1$lu;rm(d1b,d1p)

# Metadata
ab.out2 <- subset(out,resp_var=="LogAbund")
sp.out2 <- subset(out,resp_var=="Species_richness")

ab.out2$factor_levels <- droplevels(ab.out2$factor_levels)
ab.out2 <- ab.out2[which(!is.na(ab.out2$factor_levels)),]
ab.out2$land_use <- droplevels(ab.out2$land_use)
ab.out2$land_use <- ordered(ab.out2$land_use,levels=(ord))
ab.out2$upper_est[which(ab.out2$trans_est==100&ab.out2$lower_est==100&ab.out2$upper_est==100)] <- NA
ab.out2$lower_est[which(ab.out2$trans_est==100&ab.out2$lower_est==100&ab.out2$upper_est==100)] <- NA
ab.out2$Type <- "PREDICTS"

sp.out2$factor_levels <- droplevels(sp.out2$factor_levels)
sp.out2 <- sp.out2[which(!is.na(sp.out2$factor_levels)),]
sp.out2$land_use <- droplevels(sp.out2$land_use)
sp.out2$land_use <- ordered(sp.out2$land_use,levels=(ord))
sp.out2$upper_est[which(sp.out2$trans_est==100&sp.out2$lower_est==100&sp.out2$upper_est==100)] <- NA
sp.out2$lower_est[which(sp.out2$trans_est==100&sp.out2$lower_est==100&sp.out2$upper_est==100)] <- NA
sp.out2$Type <- "PREDICTS"

d2 <- rbind(ab.out2,sp.out2)
d2$Grouping <- NA
d2$model <- NULL; d2$t_or_zval <- NULL; d2$pval <- NULL; d2$use_intensity <- NULL;d2$use_intensity2 <- NULL;d2$ref_level <- NULL;d2$explan_var <- NULL;d2$factor_levels <- NULL
d <- join(d1,d2,type="full",match="all")
d$resp_var <- factor(d$resp_var);levels(d$resp_var) = c("LogAbundance","Species Richness")
d$land_use <- droplevels(d$land_use)
names(d);nrow(d)

#d$ind <- NA
#d$ind[which(d$Type=="Fieldwork")] <- paste0(d$studyID[which(d$Type=="Fieldwork")],"-",d$Grouping[which(d$Type=="Fieldwork")])

rect_lab <- c("PV","SV","PL","CL","UR")
cols<-c("#00AE00","#94BD5E", "#006B6B", "#E6E64C","#808080")
xmin = c(0.5,1.5,2.5,3.5,4.5);xmax=c(1.5,2.5,3.5,4.5,5.5)
library(ggthemr)
#ggthemr("greyscale",layout = "clean")
ggthemr("pale",layout = "clean",type = "outer",spacing = 0)
#theme_set(theme_classic(base_size=12,base_family = "sans"))
p <- ggplot(d[d$Type=="PREDICTS",],aes(x=land_use,y=trans_est,ymin=lower_est,ymax=upper_est,group=studyID)) 
# Normal
p <- p + geom_pointrange(aes(ymax=upper_est,ymin=lower_est),alpha=.8,size=.85,position=position_dodge(width=.5,height=.5) )
p <- p + geom_vline(xintercept=xmax[-5],color="black",size=.2,linetype="dotted") + geom_vline(xintercept=0)# and vline seperation
p <- p + geom_hline(y=100) + scale_size_area()
# Highlight field data
p <- p + geom_pointrange(data = d[d$Type=="Fieldwork"&d$Grouping=="Birds",],aes(ymax=upper_est,ymin=lower_est,color=studyID,fill=studyID,group=Grouping),size=1.5,position=position_jitterdodge(dodge.width=.75) )
#p <- p + geom_pointrange(data = d[d$Type=="Fieldwork"&d$Grouping=="Plants",],aes(ymax=upper_est,ymin=lower_est,color=ind,group=Grouping,fill=ind),size=1.2,position=position_jitterdodge(dodge.width = .75) )
#p <- p + geom_pointrange(data = d[which(!is.na(d$ind)),],aes(ymax=upper_est,ymin=lower_est,group=ind,color=ind,fill=ind),size=1.2,position=position_jitterdodge(dodge.width = .75) )
#p <- p + scale_colour_manual(values = cols,guide=guide_legend(title="Land use"))
p <- p + scale_fill_discrete(guide="none")
p <- p + scale_color_brewer(palette = "Set1",guide=guide_legend(title="Independent data",nrow = 2)) # Remove legend from jitterdodge
p <- p + scale_x_discrete(labels=rect_lab) 
p <- p + facet_grid(~resp_var)
p <- p + coord_cartesian(ylim=c(0,400)) 
p <- p + scale_y_continuous(expand=c(1,1),breaks=pretty_breaks(),limits=c(0,400)) 
p <- p + labs(x="",y="Percentage of Baseline Diversity Metric")
p <- p + theme(axis.ticks = element_blank(),
                axis.text.y=element_text(size=10),
                axis.title.y = element_text(size=16,vjust=0.5,face="bold"),
                strip.text.x = element_text(size = 20),
                legend.text=element_text(size=14),                
                legend.position="bottom")
#p <- p + scale_shape_manual(breaks = levels(factor(d$Grouping[d$Type=="Fieldwork"])),guide=guide_legend(title="Fieldwork Datasets"),
#                   values = c(17, 19))
p
ggsave(paste("OutputComparison/Consistency_facet.png"),plot=p,width = 12,height=9,dpi=400)
#FOR THESIS
ggsave("../000Thesis_writeup/WriteUpLatex/gfx/SIFigure2.pdf",plot=p,scale=1.1,dpi=400)

# Multiplot
sites.pred$LogAbund <-sites.pred$logabund
sites.pred$ChaoRint <-sites.pred$ChaoR
sites.pred$LandUse <- sites.pred$PREDICTS.LU
sites.pred$UseIntensity <- sites.pred$Use_intensity
sites.pred$UI <- as.factor(paste(sites.pred$LandUse,sites.pred$UseIntensity))

#sub.sites.pred <- subset(sites.pred,PREDICTS.LU!=c("Urban"))
sub.sites.pred <- subset(sites.pred,PREDICTS.LU!=c("Pasture"))
sub.sites.pred$PREDICTS.LU <-droplevels(sub.sites.pred$PREDICTS.LU)

fit <- glm(logabund~PREDICTS.LU,data=sub.sites,family=gaussian)
fit <- glm(spec~CroplandSplit,data=sub.sites,family=poisson(link=log))
pfit <- glmer(LogSimpson~LandUse+(1|SS),data=sub.sites.pred,family=gaussian)
pfit <- glmer(LogAbund~LandUse+(1|SS),data=sub.sites.pred,family=gaussian)
pfit <- glmer(Species_richness~CroplandSplit+(1|SS),data=sub.sites.pred,family=poisson(link=log))

up <- stepAIC(fit,direction="both",trace=T)#up = fit
summary(up)

plot(Effect(focal.predictors="CroplandSplit",mod=up),ylab="Species Richness",xlab="LandUse",main="Fieldwork glm")
plot(Effect(focal.predictors="CroplandSplit",mod=pfit),ylab="Species Richness",xlab="LandUse",main="PREDICTS glmer")

par(pty="s",las=1, bty="l",pch=19, cex=1.0,cex.lab=1,cex.axis=0.7, tcl=-0.2, mar=c(7,1,3,0)+0.1, mgp=c(2.5,0.5,0))
boxplot(trans_est~land_use,data=dfc,las=2);abline(h=0,lty=2,lwd=.5)
par(par.ori)

#### PREDICTS Full model comparison ####

#cbind(sum(residuals(pfit_lu_simple,type="pearson")^2),df.residual(pfit_lu_simple))
#cbind(sum(residuals(pfit_lu_simple_rich,type="pearson")^2),df.residual(pfit_lu_simple_rich))
# Sum of squared Pearson residuals is way less than the residual degrees of freedom,
# thus not overdispersed

library(marfunky)
standardPackages("Analysis")
sem <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

sub.sites <- sites#subset(sites,PREDICTS.LU!="Urban")
sub.sites$LogAbund <- sub.sites$logabund
sub.sites$LogSimpson <- sub.sites$logsimp
sub.sites$Species_richness <- sub.sites$spec

sub.sites.pred <- subset(sites.pred,PREDICTS.LU!=c("Pasture"))
sub.sites.pred$PREDICTS.LU <-droplevels(sub.sites.pred$PREDICTS.LU)
sub.sites.pred_b <- subset(sub.sites.pred,Grouping==c("Birds")) # Get only birds
sub.sites.pred_b$PREDICTS.LU <-droplevels(sub.sites.pred_b$PREDICTS.LU)
sub.sites.pred_p <- subset(sub.sites.pred,Grouping==c("Plants")) # Get only plants
sub.sites.pred_p$PREDICTS.LU <-droplevels(sub.sites.pred_p$PREDICTS.LU)

#(table(sites.pred$Grouping))
makeRatios <- function(model,pred="PREDICTS.LU",exp=T,ci=F,boot=100){
  fix <- fixef(model)
  se.fix <- se.fixef(model)
#   if(!("(Intercept)" %in% names(fix))){
#     se.fix <- c(se.fixef(model)[1],se.fixef(model)[2:5]-se.fixef(model)[1])
#   }

  if("(Intercept)" %in% names(fix)){
    df <- data.frame(LandUse=c("Primary Vegetation",sub(pred,replacement = "",names(fix)[-1])),
                     eff=NA,se.low=NA,se.high=NA,Transect="PREDICTS")    
    if(exp){
      df$eff <-(exp(fix))*100
      df$se.high <-(exp(fix+se.fix))*100
      df$se.low <-(exp(fix-se.fix))*100
    } else {
      df$eff <-((fix))*100
      df$se.high <-((fix+se.fix))*100
      df$se.low <-((fix-se.fix))*100
    }
    # Reset the intercept to Zero
    df[1,2:4] <- c(100,100,100)
    df$LandUse <- ordered(df$LandUse,levels=(df$LandUse))
    # Using Confidence intervals or se
    if(ci) {
      df$ci.low=NA;df$ci.high=NA
      cf <- confint.merMod(model,quiet = T,level = .95,method="boot",nsim = boot,
                                  .progress="txt",PBargs=list(style=3))
      
      df$ci.high <- (exp(cf[3,2]+cf[names(fix),1]))
      df$ci.low <- (exp(cf[3,1]-cf[names(fix),2]))
      df[c(1,5),] <- 100
    }
    } else {
    

    df <- data.frame(LandUse=sub(pred,replacement = "",names(fix)),
                     eff=NA,se.low=NA,se.high=NA,Transect="PREDICTS")
    df$LandUse <- ordered(df$LandUse,levels=(df$LandUse))
    # Using Confidence intervals or se
    if(ci){  cf <- confint.merMod(model,quiet = T,method="boot",nsim = boot,
                                  .progress="txt",PBargs=list(style=3))
    } 
    # Get Standard error
    if(exp==T) {

    se.high <- exp(fix + se.fix )
    se.low <- exp(fix - se.fix )
    fix <- exp(fix)
    } else{
      se.high <- (fix + se.fix )
      se.low <- (fix - se.fix )
    }
    if(ci) {df$ci.low=NA;df$ci.high=NA}
    # Correct errors for baseline percentage
    for(i in 1:length(fix)){    
      df$eff[i] <- (fix[i] / fix[1]) * 100
      df$se.low[i] <- (se.low[names(fix)[i]] / se.low[names(fix)[1]]) * 100
      df$se.high[i] <- (se.high[names(fix)[i]]  / se.high[names(fix)[1]]) * 100
      if(ci){
        df$ci.low[i] <- (cf[names(fix)[i],1] / cf[names(fix)[1],1]) * 100
        df$ci.high[i] <- (cf[names(fix)[i],2] / cf[names(fix)[1],2]) * 100
      }
    } 
  }
  return(df)
}

rc <- compare_randoms(sub.sites.pred,"logabund","gaussian",fixedFactors=c("PREDICTS.LU"),siteRandom=F)
pfit_lu_simple <- glmer(as.formula(paste("logabund ~  PREDICTS.LU +",rc$best.random)),data=sub.sites.pred,family=gaussian,REML=F)
rc <- compare_randoms(sub.sites.pred,"Species_richness","poisson",fixedFactors="PREDICTS.LU",siteRandom=T)
pfit_lu_simple_rich <- glmer(as.formula(paste("Species_richness ~ PREDICTS.LU + ",rc$best.random)),data=sub.sites.pred,family=poisson(link=log),REML=F)
dfc <- makeRatios(pfit_lu_simple,"PREDICTS.LU",T,F);dfc$data <- "All";dfc$type <- "log(Abundance)";dfc$Grouping <- NA
dfc_rich <- makeRatios(pfit_lu_simple_rich,"PREDICTS.LU",T,F);dfc_rich$data <- "All";dfc_rich$type <- "Species Richness";dfc_rich$Grouping <- NA
# Birds only
rc <- compare_randoms(sub.sites.pred_b,"logabund","gaussian",fixedFactors="PREDICTS.LU",siteRandom=F)
pfit_lu_simple <- glmer(as.formula(paste("logabund ~  PREDICTS.LU +",rc$best.random)),data=sub.sites.pred_b,family=gaussian,REML=T)
rc <- compare_randoms(sub.sites.pred,"Species_richness","poisson",fixedFactors="PREDICTS.LU",siteRandom=F)
pfit_lu_simple_rich <- glmer(as.formula(paste("Species_richness ~ PREDICTS.LU + ",rc$best.random)),data=sub.sites.pred_b,family=poisson(link=log),REML=T)
dfc_b <- makeRatios(pfit_lu_simple,"PREDICTS.LU",T,F);dfc_b$data <- "Birds";dfc_b$type <- "log(Abundance)";dfc_b$Grouping <- NA
dfc_rich_b <- makeRatios(pfit_lu_simple_rich,"PREDICTS.LU",T,F);dfc_rich_b$data <- "Birds";dfc_rich_b$type <- "Species Richness";dfc_rich_b$Grouping <- NA
# Plants only
#pfit_lu_simple <- glmer(logabund ~ PREDICTS.LU + (1|SS),data=sub.sites.pred_b,family=gaussian,REML=T)
#pfit_lu_simple_rich <- glmer(Species_richness ~ PREDICTS.LU + (1|SS),data=sub.sites.pred_b,family=poisson(link=log),REML=T)
#dfc_p <- makeRatios(pfit_lu_simple,"PREDICTS.LU",T,F);dfc_p$data <- "Plants";dfc_p$type <- "log(Abundance)";dfc_p$Grouping <- NA
#dfc_rich_p <- makeRatios(pfit_lu_simple_rich,"PREDICTS.LU",T,F);dfc_rich_p$data <- "Plants";dfc_rich_p$type <- "Species Richness";dfc_rich_p$Grouping <- NA

# Add an NA row for non-existing urban for now
#dfc_b <- rbind(dfc_b,data.frame(LandUse="Urban",eff=NA,se.low=NA,se.high=NA,Transect="PREDICTS",data="Birds",type="log(Abundance)",Grouping=NA))
#dfc_rich_b <- rbind(dfc_rich_b,data.frame(LandUse="Urban",eff=NA,se.low=NA,se.high=NA,Transect="PREDICTS",data="Birds",type="Species Richness",Grouping=NA))

# Site Data
# Build indices for each Grouping
sA <- subset(sub.sites,Grouping=="Aves")
sP <- subset(sub.sites,Grouping=="Plantae")

sitemakeRatios <- function(sub.sites,pred,exp=T,Grouping="Birds") {
  fw_ab <- summarySEwithin(sub.sites,measurevar="val",withinvars=c("Transect",pred))
  fw_ab$ind <- match(fw_ab[,2],ord)
  fw_ab$se.high <- fw_ab$val + fw_ab$se  
  fw_ab$se.low <- fw_ab$val - fw_ab$se
  if(exp){
    fw_ab$val <- exp(fw_ab$val)
    fw_ab$se.high <- exp(fw_ab$se.high)
    fw_ab$se.low <- exp(fw_ab$se.low)
  }
  # Do for Taita
  if ("Primary Vegetation" %in% unique(subset(fw_ab,Transect=="Taita")[,pred])){
    LUrefT <- 1
  } else {LUrefT <- 2}
  t <- fw_ab$val[which(fw_ab$ind==LUrefT&fw_ab$Transect=="Taita")]
  if ("Primary Vegetation" %in% unique(subset(fw_ab,Transect=="Kilimanjaro")[,pred])){
    LUrefK <- 1
  } else {LUrefK <- 2} 
  k <- fw_ab$val[which(fw_ab$ind==LUrefK&fw_ab$Transect=="Kilimanjaro")]
  
  fw_ab$nratio[which(fw_ab$Transect=="Taita")] <- (fw_ab$val[which(fw_ab$Transect=="Taita")]/t )*100
  fw_ab$nratio[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$val[which(fw_ab$Transect=="Kilimanjaro")]/k)*100
  

  fw_ab$se.low.n[which(fw_ab$Transect=="Taita")] <- (fw_ab$se.low[which(fw_ab$Transect=="Taita")]/fw_ab$se.low[which(fw_ab$ind==LUrefT&fw_ab$Transect=="Taita")] )*100
  fw_ab$se.high.n[which(fw_ab$Transect=="Taita")] <- (fw_ab$se.high[which(fw_ab$Transect=="Taita")]/fw_ab$se.high[which(fw_ab$ind==LUrefT&fw_ab$Transect=="Taita")])*100
  fw_ab$se.low.n[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$se.low[which(fw_ab$Transect=="Kilimanjaro")]/fw_ab$se.low[which(fw_ab$ind==LUrefK&fw_ab$Transect=="Kilimanjaro")] )*100
  fw_ab$se.high.n[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$se.high[which(fw_ab$Transect=="Kilimanjaro")]/fw_ab$se.high[which(fw_ab$ind==LUrefK&fw_ab$Transect=="Kilimanjaro")])*100
  
  of <- data.frame(LandUse=fw_ab[,2],eff=fw_ab$nratio,se.low=fw_ab$se.low.n,se.high=fw_ab$se.high.n,Transect=fw_ab$Transect,data="Field",Grouping=Grouping)
  of$LandUse <- droplevels(ordered(of$LandUse,levels=(ord)))
  return(of)
}

detachPackage("dplyr")
sA$val <- sA$logabund;sP$val <- sP$logabund
ofA <- sitemakeRatios(sA,"PREDICTS.LU",exp=T,"Birds")
ofP <- sitemakeRatios(sP,"PREDICTS.LU",exp=T,"Plants")
ofA$type <- "log(Abundance)";ofP$type <- "log(Abundance)"

sA$val <- sA$spec;sP$val <- sP$spec
ofA_r <- sitemakeRatios(sA,"PREDICTS.LU",exp=F,"Birds")
ofP_r <- sitemakeRatios(sP,"PREDICTS.LU",exp=F,"Plants")
ofA_r$type <- "Species Richness";ofP_r$type <- "Species Richness"

full <- rbind(ofA,ofA_r,dfc,dfc_rich)#,ofP,ofP_r)
#full <- rbind(ofA,ofA_r,dfc_b,dfc_rich_b) # Birds only in PREDICTS
full$LandUse <- ordered(full$LandUse,ord)
full$Transect <- ordered(full$Transect,levels=c("Kilimanjaro","PREDICTS","Taita"))
#full$NewGroup <- paste0(full$Transect,ifelse(is.na(full$Grouping),"",paste0("-",full$Grouping)))
#full$NewGroup <- ordered(full$NewGroup,levels=c("Kilimanjaro-Plants","Kilimanjaro-Birds","PREDICTS","Taita-Birds","Taita-Plants"))
#full$Grouping2 <- ifelse(is.na(full$Grouping),"Birds",as.character(full$Grouping))
#f <- full[which(full$Transect=="PREDICTS"),];f$Grouping2 <- "Plants"; full <- rbind(full,f);rm(f)

### Create Z- scores ###
library(dplyr)
library(tidyr)
full2 <- full %>% filter(type == "log(Abundance)") %>% 
  mutate(se = abs(se.high - eff)) %>%  #recreate se
  dplyr::select(LandUse, eff,se, data) %>% # deselect
  group_by(LandUse,data) %>% 
  summarize(eff = mean(eff),
            se = mean(se)) %>% ungroup() # summarize
d1 = full2 %>% select(-se) %>% spread(data,value = eff) %>% rename(Field.eff = Field, All.eff = All) 
d2 = full2 %>% select(-eff) %>% spread(data,value = se) %>% rename(Field.se = Field, All.se = All) 
d_a <- left_join(d1,d2)  %>% 
  #Calculate Z score
  mutate(Z = (Field.eff - All.eff) / sqrt((All.se * All.se) + (Field.se * Field.se)),
         type = "log(Abundance)") %>% dplyr::select(LandUse,type,Z)

full2 <- full %>% filter(type == "Species Richness") %>% 
  mutate(se = abs(se.high - eff)) %>%  #recreate se
  dplyr::select(LandUse, eff,se, data) %>% # deselect
  group_by(LandUse,data) %>% 
  summarize(eff = mean(eff),
            se = mean(se)) %>% ungroup()  # summarize
d1 = full2 %>% select(-se) %>% spread(data,value = eff) %>% rename(Field.eff = Field, All.eff = All) 
d2 = full2 %>% select(-eff) %>% spread(data,value = se) %>% rename(Field.se = Field, All.se = All) 
d_s <- left_join(d1,d2)  %>% 
  #Calculate Z score
  mutate(Z = (Field.eff - All.eff) / sqrt((All.se * All.se) + (Field.se * Field.se)),  
         type = "Species Richness") %>% dplyr::select(LandUse,type,Z)
dsa <- rbind(d_a,d_s) %>% mutate(Transect = NA)
dsa <- dsa[-c(1,6),]
# --------------- #
dsa %>% group_by(type) %>% summarise(Zm = median(abs(Z),na.rm=T),
                                     Zmin = min(abs(Z),na.rm=T),
                                     Zmax = max(abs(Z),na.rm=T))

# Pub plot
rect_lab <- c("PV","SV","PL","CL","UR")
xmin = c(0.5,1.5,2.5,3.5,4.5);xmax=c(1.5,2.5,3.5,4.5,5.5)
theme_set(theme_classic(base_size=12,base_family = "sans"))
g <- ggplot(full,aes(x=LandUse,y=eff,shape=Transect,group=Transect)) + theme_classic(base_size = 16)
g <- g + geom_hline(yintercept=100,color="lightblue", lty=2,size=.5)
g <- g + geom_vline(xintercept=xmax[-5],size=.2,linetype="dotted") + geom_vline(xintercept=0)# and vline seperation
#g <- g + geom_vline(xintercept=xmax[-5],size=.2,color="white",linetype="dotted") + geom_vline(xintercept=0)# and vline seperation
#g <- g + geom_pointrange(aes(ymax=se.high,ymin=se.low),size=.9,position=position_dodge(width=1) )
g <- g + geom_pointrange(aes(ymax=se.high,ymin=se.low),size=.9,position=position_dodge(width=1) )
g <- g + facet_wrap(~type,scales = "free_y",as.table = F,drop = T) #+ theme(strip.text.x = element_text(size=16), strip.background = element_rect(colour="black"))
g <- g + scale_y_continuous(expand=c(0,1),breaks=pretty_breaks()) + coord_cartesian()
#g <- g + theme(axis.text.y=element_text(face = c(rep("plain",4),"bold",rep("plain",10)) )) # Special Eyecandy of 100% 
g <- g + scale_x_discrete(labels=rect_lab)
# Insert z-score
g <- g + geom_label(data=dsa,aes(x = LandUse, y = c(rep(177,4),rep(110,4)),label = paste("Z = ",round(Z,2) )))
g <- g + labs(y="Percentage of Baseline",x="")
#g <- g + scale_shape_manual(values=c(15,16,17),labels = c("Kilimanjaro","Broad-scale","Taita"),
#                            guide = guide_legend(title = "Dataset:",direction = "horizontal", title.position = "left",
#                                                 label.position="right", label.hjust = 0.5, label.vjust = 0.5,nrow=1)) 
#g <- g + theme(axis.text.x = element_text(size = rel(1.2), colour = 'black')) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
#g <- g + theme(axis.title.x = element_text(size =rel(1.3), colour = 'black')) + theme(axis.title.y = element_text(size = 20, colour = 'black'))
g <- g + theme(legend.position="bottom",axis.title.x=element_blank(),axis.ticks.x=element_blank())
g
ggsave("../MS_Submission/Round2/Appendix/Figure2.png",plot=g,units = "mm",dpi=400)
#ggsave("../Conferences/Conservation Symposium/Comparison.png",plot=g,scale = 1.1,units = "mm",dpi=400)
# FOR THESIS
ggsave("../000Thesis_writeup/WriteUpLatex/gfx/Figure2.pdf",plot=g,scale = 1.1,units = "mm",dpi=400)


# Transparent
t <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "transparent",colour = NA),
           axis.text = element_text(colour="black"),axis.text.x =element_text(colour="black"),axis.text.y =element_text(colour="black"),axis.ticks = element_line(colour = "black"), axis.line = element_line(colour="black"),axis.title = element_text(colour="black"),
           legend.key = element_rect(fill=NA,colour = NA),legend.background = element_rect(fill = NA,colour=NA),legend.text = element_text(colour="black"),
           title = element_text(colour="black",size=18),
           strip.background = element_rect(fill=NA),
           panel.grid.major = element_blank(),panel.grid.minor = element_blank()
) 
ggsave("../001Thesis_Defense/assets/img/Figure2_trans.png",plot=g+t,scale = 1.1,units = "mm",dpi=400,bg="transparent")

# PREDICTS bird only comparison
library(dplyr);detachPackage("MODIS");detachPackage("rasterVis");detachPackage("raster")
full_pred <- rbind(dfc,dfc_rich,dfc_b,dfc_rich_b)
full_pred$LandUse <- ordered(full_pred$LandUse,ord)
# Points of special Interest
poi <- rbind(ofA,ofA_r) %>% dplyr::filter(LandUse=="Cropland")

rect_lab <- c("PV","SV","PL","CL","UR")
xmin = c(0.5,1.5,2.5,3.5,4.5);xmax=c(1.5,2.5,3.5,4.5,5.5)
ggthemr("greyscale",layout = "clean",type="outer",spacing=1)
g <- ggplot(full_pred,aes(x=LandUse,y=eff,shape=data,group=data))
g <- g + geom_hline(yintercept=100,color="lightblue", lty=2,size=.5)
g <- g + geom_vline(xintercept=xmax[-5],size=.2,color="white",linetype="dotted") + geom_vline(xintercept=0)# and vline seperation
g <- g + geom_pointrange(aes(ymax=se.high,ymin=se.low),size=.9,position=position_dodge(width = 0.5) )
g <- g + facet_wrap(~type,scales = "free_y",nrow = 1,as.table = F,drop = T)
g <- g + scale_y_continuous(expand=c(0,1),breaks=pretty_breaks()) + coord_cartesian()
g <- g + scale_x_discrete(labels=rect_lab)
g <- g + labs(y="Percentage of Baseline",x="")
g <- g + scale_shape_manual(values=c(21,25,19),labels = c("All taxonomic groups","Only birds","Cropland field transects"),
                            guide = guide_legend(title = "Dataset:",direction = "horizontal", title.position = "left",
                                                 label.position="right", label.hjust = 0.5, label.vjust = 0.5,nrow=1)) 
g <- g + theme(legend.position="bottom",axis.title.x=element_blank(),axis.ticks.x=element_blank())
# Add POI
g <- g + geom_pointrange(data=poi,mapping = aes(ymax=se.high,ymin=se.low,x=LandUse,y=eff,group=Transect),size=1.2,color="black",position=position_dodge(width=1),show_guide = FALSE)
g
ggsave("../000Thesis_writeup/PaperDraft/Appendix/SI Figure3.png",plot=g,scale = 1.1,units = "mm",dpi=400)

# ------------------- #
rect_lab <- c("PV","SV","PL","CL","UR")
xmin = c(0.5,1.5,2.5,3.5,4.5);xmax=c(1.5,2.5,3.5,4.5,5.5)
library(ggthemr)
#ggthemr("earth",layout = "scientific",type = "outer",spacing = 0)
#ggthemr("greyscale",layout = "clean")
ggthemr("fresh",layout="clean",type="outer",spacing=0)
#theme_set(theme_classic(base_size=12,base_family = "sans"))
g <- ggplot(full,aes(x=LandUse,y=eff,shape=Transect,group=Transect,color=Transect))
g <- g + geom_hline(yintercept=100,color="lightblue", lty=2,size=.5)
#g <- g + geom_vline(xintercept=xmax[-5],size=.2,linetype="dotted") + geom_vline(xintercept=0)# and vline seperation
g <- g + geom_vline(xintercept=xmax[-5],size=.2,color="white",linetype="dotted") + geom_vline(xintercept=0)# and vline seperation
#g <- g + geom_pointrange(aes(ymax=se.high,ymin=se.low),size=.9,position=position_dodge(width=1) )
g <- g + geom_pointrange(aes(ymax=se.high,ymin=se.low),size=.9,position=position_dodge(width=1) )

g <- g + facet_wrap(~type,scales = "free_y",as.table = F,drop = T) #+ theme(strip.text.x = element_text(size=16), strip.background = element_rect(colour="black"))
g <- g + scale_y_continuous(expand=c(0,1),breaks=pretty_breaks()) + coord_cartesian()
#g <- g + theme(axis.text.y=element_text(face = c(rep("plain",4),"bold",rep("plain",10)) )) # Special Eyecandy of 100% 
g <- g + scale_x_discrete(labels=rect_lab)
g <- g + labs(y="Percentage of Baseline",x="")
#g <- g + scale_shape_manual(values=c(15,16,17),labels = c("Kilimanjaro","Broad-scale","Taita"),
#                            guide = guide_legend(title = "Dataset:",direction = "horizontal", title.position = "left",
#                                                 label.position="right", label.hjust = 0.5, label.vjust = 0.5,nrow=1)) 
#g <- g + theme(axis.text.x = element_text(size = rel(1.2), colour = 'black')) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
#g <- g + theme(axis.title.x = element_text(size =rel(1.3), colour = 'black')) + theme(axis.title.y = element_text(size = 20, colour = 'black'))
g <- g + theme(legend.position="bottom",axis.title.x=element_blank(),axis.ticks.x=element_blank())
g
ggsave("OutputComparison/Comparison_pubready.png",plot=g,scale = 1.5,units = "mm",dpi=400)




# ---------------------- #
rect_colours<-c("#00AE00","#94BD5E", "#006B6B", "#E6E64C","#808080")
rect_lab <- c("PV","SV","PL","CL","UR")
xmin = c(0.5,1.5,2.5,3.5,4.5);xmax=c(1.5,2.5,3.5,4.5,5.5)
rem <- theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank()) # Remove legend and x-axis stuff
#rem_la <- theme(axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank()) # Remove left axis and labels as well 

# With rectangulars
#d <- rbind(dfc,of)
full$LandUse <- droplevels(full$LandUse)
levels(full$LandUse) <- rect_lab
theme_set(theme_classic(base_size=24,base_family = "sans"))
full1 <- subset(full,type=="log(Abundance)")
g <- ggplot(full1,aes(x=LandUse,y=eff,color=Transect,group=Transect))
g <- g + geom_pointrange(aes(ymax=se.high,ymin=se.low),size=1,position=position_dodge(width=.5) )
# Do the same with grids
g <- g + annotate("rect", xmin = xmin,xmax = xmax,ymin = -Inf, ymax = Inf,fill = rect_colours,alpha = rep(.5,5)) 
g <- g + geom_hline(yintercept=100, lty=2,size=.5)
g <- g + geom_pointrange(aes(ymax=se.high,ymin=se.low),size=1,position=position_dodge(width=.5) )
g <- g + labs(y ='Percentage of Baseline', x = "", title=unique(full1$type))
g <- g + scale_y_continuous(expand=c(0,0),limits=c(50,185),breaks = pretty_breaks(10))
g <- g + scale_color_manual(values=c("red","black","blue"),labels = c("Kilimanjaro","PREDICTS","Taita Hills"),
                            guide = guide_legend(direction = "horizontal", title.position = "left",title="",
                                                 label.position="right", label.hjust = 0.5, label.vjust = 0.5)) 
g <- g + theme(axis.text.x = element_text(size = rel(1.2), colour = 'black')) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
g <- g + theme(axis.title.x = element_text(size =rel(1.3), colour = 'black')) + theme(axis.title.y = element_text(size = 20, colour = 'black'))
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g <- g + theme(legend.position="none") + theme(title = element_text(size=22)) + rem

full2 <- subset(full,type=="Species Richness")
levels(full2$LandUse) <- c("Primary\n vegetation","Secondary\n vegetation","Plantation\n forest","Cropland","Urban")
full2$Transect <- factor(full2$Transect, levels = c("Kilimanjaro","Taita","PREDICTS"))
g2 <- ggplot(full2,aes(x=LandUse,y=eff,color=Transect,shape=data,group=Transect))
g2 <- g2 + geom_pointrange(aes(ymax=se.high,ymin=se.low),size=2,position=position_dodge(width=.5) )
# Do the same with grids
g2 <- g2 + annotate("rect", xmin = xmin,xmax = xmax,ymin = -Inf, ymax = Inf,fill = rect_colours,alpha = rep(.5,5)) 
g2 <- g2 + geom_hline(yintercept=100, lty=2,size=.5)
g2 <- g2 + geom_pointrange(aes(ymax=se.high,ymin=se.low),size=2,position=position_dodge(width=.5) )
g2 <- g2 + labs(y =paste0('Species richness \n relative to primary vegetation ',"(",expression("%"),")"), x = "", title="")
g2 <- g2 + scale_y_continuous(expand=c(0,0),limits=c(50,120),breaks = pretty_breaks())
g2 <- g2 + scale_shape_manual(values=c(16,15))
g2 <- g2 + scale_color_manual(values=c("red","blue","black"),labels = c("Kilimanjaro","PREDICTS","Taita Hills"),
                            guide = guide_legend(direction = "horizontal", title.position = "left",title="",
                                                 label.position="right", label.hjust = 0.5, label.vjust = 0.5)) 
g2 <- g2 + theme(axis.text.x = element_text(size = rel(1.2), colour = 'black')) +theme(axis.text.y = element_text(size = 24, colour = 'black'))
g2 <- g2 + theme(axis.title.x = element_text(size =rel(1.3), colour = 'black')) + theme(axis.title.y = element_text(size = 26, colour = 'black'))
g2 <- g2 + theme(axis.text.x = element_text(angle = 0, size=24),axis.ticks.x = element_blank())
g2 <- g2 + theme(legend.position="none") + theme(title = element_text(size=24)) 

library(cowplot)
t <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "transparent",colour = NA),
           axis.text = element_text(colour="black"),axis.text.y =element_text(colour="black"), axis.line = element_line(colour="black"),axis.title = element_text(colour="black"),
           legend.key = element_rect(fill=NA,colour = NA),legend.background = element_rect(fill = NA,colour=NA),legend.text = element_text(color="black"),
           title = element_text(colour="black"),
           strip.background = element_rect(fill=NA,color=NA), strip.text = element_text(color="black"),
           panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
           text = element_text(color="black",size=24)
) 
gg <- plot_grid(g+t,g2+t)

ggplot2::ggsave("Figure2_SR_trans.png",plot=g2+t,dpi=400,scale = 1.1,units = "mm",bg="transparent")


"Species Richness"
# -------- #
r.squaredGLMM(pfit_lu_simple)
r.squaredGLMM(pfit_lu_simple_rich)

# Stats #
full <- rbind(ofA,ofA_r,dfc,dfc_rich)#,ofP,ofP_r)
full_b <- rbind(ofA,ofA_r,dfc_b,dfc_rich_b)
#full_p <- rbind(ofP,ofP_r,dfc_p,dfc_rich_p)

# Compute goodness of fit --> rsquared
a <- as.data.frame(acast(full, type+LandUse~Transect, value.var="eff",fun.aggregate = mean))
a2 <- as.data.frame(acast(full, type+LandUse~data, value.var="eff",fun.aggregate = mean))
b <- as.data.frame(acast(full_b, type+LandUse~Transect+data, value.var="eff",fun.aggregate = mean))
b2 <- as.data.frame(acast(full_b, type+LandUse~data, value.var="eff",fun.aggregate = mean))
#p <- as.data.frame(acast(full_p, type+LandUse~Transect+data, value.var="eff",fun.aggregate = mean))
#p2 <- as.data.frame(acast(full_p, type+LandUse~data, value.var="eff",fun.aggregate = mean))

# manually computation of chi-square = (Observed - expected)² / sqrt(expected)
a$k_chisquare <- (a$Kilimanjaro-a$PREDICTS)^2 / sqrt(a$PREDICTS) 
a$t_chisquare <- (a$Taita-a$PREDICTS)^2 / sqrt(a$PREDICTS)

# ---- # 
# Publication

# Manually calculate the adjusted r². 
# 1 - [Sum(i=1 to n){wi (yi - fi)2}] /[Sum(i=1 to n){wi (yi - yav)2}]
# Or use lm? Yes, identical to squared pearson r with one predictor

# Overall
summary(lm(Field~All,data = a2,subset=1:5)) # abu
summary(lm(Field~All,data = a2,subset=6:10)) # sr

# Without crop
summary(lm(Field~All,data = a2,subset=c(1:3,5))) # abu
summary(lm(Field~All,data = a2,subset=c(6:8,10))) # sr

# Transect
summary(lm(Taita~PREDICTS,data = a,subset=1:5)) # abu
summary(lm(Taita~PREDICTS,data = a,subset=6:10)) # sr
summary(lm(Kilimanjaro~PREDICTS,data = a,subset=1:5)) #abu
summary(lm(Kilimanjaro~PREDICTS,data = a,subset=6:10)) #sr

# Individual birds
summary(lm(Taita_Field~PREDICTS_Birds,data = b,subset=1:5)) # abu
summary(lm(Taita_Field~PREDICTS_Birds,data = b,subset=6:10)) # sr
summary(lm(Kilimanjaro_Field~PREDICTS_Birds,data = b,subset=1:5)) #abu
summary(lm(Kilimanjaro_Field~PREDICTS_Birds,data = b,subset=6:10)) #sr

#summary(lm(p$Taita_Field~p$PREDICTS_Plants))
#summary(lm(p$Kilimanjaro_Field~p$PREDICTS_Plants))

library(car)
# All abundance and spec
cor.test(a2$Field[1:5],a2$All[1:5],method = "pearson",exact = T)$estimate^2
cor.test(a2$Field[6:10],a2$All[6:10],method = "pearson",exact = T)

cor.test(a$Kilimanjaro,a$PREDICTS,method = "pearson",exact = T)
cor.test(a$Taita,a$PREDICTS,method = "pearson",exact = T)

cor.test(b2$Field,b2$Birds,method = "pearson",exact = T)
#cor.test(p2$Field,p2$Plants,method = "pearson",exact = T)


#### Investigate Species specific patterns in PREDICTS model ####
sfit <- glmer(Species_richness~PREDICTS.LU + (0+Grouping|PREDICTS.LU) + (1|SS),data=sites.pred,family = poisson)
afit <- glmer(logabund~PREDICTS.LU + (0+Grouping|PREDICTS.LU) + (1|SS),data=sites.pred,family = gaussian)

a = ranef(afit)$PREDICTS.LU # variation of taxa within LandUse
a$LU <- rownames(a)
m <- melt(a,id.vars = "LU")
m$variable <- str_replace_all(m$variable,pattern = "Grouping",replacement = "")
# order
ord <- c("Primary Vegetation", "Secondary Vegetation","Plantation forest","Pasture", "Cropland","Urban")
m$LU <- ordered(m$LU,ord)
gr <- rev(c("Mammals","Birds","Amphibia","Reptiles","Invertebrates","Plants","Other"))
m$variable <- ordered(m$variable,gr)

theme_set(theme_classic(base_size=12,base_family = "sans"))
g <- ggplot(m, aes(LU, variable, fill = value))
g <- g + geom_tile()
#g <- g + scale_fill_gradientn("Random Slope Effect",colours=heat.colors(5))
g <- g + scale_fill_gradient2("Effect",low = "red",mid = "white",high = "green",midpoint = 0)
g <- g + ggtitle("Taxon specific random slope") + ylab("") + xlab("")
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g

# Plot facets per Group and random slope effect on Landuse?
#### Communtiy overlap ####
# Compare Community overlap between higher diversity area(primary) and cropland 
# In predicts data use the intercept of ordered classes
# first check roquefort for sorensen
# present as average community overlap between high-diversity and cropland habiat to show
# quality of homegardens -> higher similarity to forest

#Sorensen = bray on binary
modCompDiss <- function (data, metric,crop=T) {
  if (metric == "SorAbd") {
    data <- data[data$Diversity_metric_type == "Abundance", 
                 ]
  }
  if (metric == "SorCorr") {
    data <- data[((data$Diversity_metric == "Abundance") & 
                    (data$Diversity_metric_unit == "individuals")), ]
  }
  data$LandUse<-paste(data$Predominant_habitat)
  data$LandUse[which(data$LandUse=="Primary forest")]<-"Primary Vegetation"
  data$LandUse[which(data$LandUse=="Primary non-forest")]<-"Primary Vegetation"
  data$LandUse[which(data$LandUse=="Young secondary vegetation")]<-"Secondary Vegetation"
  data$LandUse[which(data$LandUse=="Intermediate secondary vegetation")]<-"Secondary Vegetation"
  data$LandUse[which(data$LandUse=="Mature secondary vegetation")]<-"Secondary Vegetation"
  data$LandUse[which(data$LandUse=="Secondary vegetation (indeterminate age)")]<-"Secondary Vegetation"
  data$LandUse[which(data$LandUse=="Secondary non-forest")]<-"Secondary Vegetation"
  data$LandUse[which(data$LandUse=="Cannot decide")]<-NA
  data$LandUse<-factor(data$LandUse)
  
  #data$LandUse <- factor( data$LandUse, c("Primary Vegetation", "Secondary Vegetation","Plantation forest","Pasture", "Cropland","Urban"))
  #data$LandUse<-relevel(data$LandUse,ref="Primary Vegetation")
  data <- subset(data, select = c("SS", "SSBS", "Measurement", 
                                  "Taxon_name_entered", "LandUse"))
  if(crop){
    data <- subset(data,LandUse!="Pasture")
    data$LandUse <- droplevels(data$LandUse)
  }
  data$LandUse <- ordered( data$LandUse, c("Primary Vegetation", "Secondary Vegetation","Plantation forest","Pasture", "Cropland","Urban"))  
  data <- na.omit(data)
  if (metric == "SorCorr") {
    study.all.int.meas <- tapply(data$Measurement, data$SS, 
                                 function(m) all(floor(m) == m))
    int.meas <- study.all.int.meas[match(data$SS, names(study.all.int.meas))]
    data <- data[int.meas, ]
  }
  all.lu <- unique(paste(data$LandUse))
  all.results <- matrix(nrow = length(all.lu), ncol = length(all.lu))
  sum.matrix <- matrix(0, nrow = length(all.lu), ncol = length(all.lu))
  count.matrix <- matrix(0, nrow = length(all.lu), ncol = length(all.lu))
  all.results <- data.frame(all.results)
  names(all.results) <- all.lu
  row.names(all.results) <- all.lu
  all.studies.count <- 0
  used.studies <- 0
  for (st in unique(data$SS)) {
    all.studies.count <- all.studies.count + 1
    cat(paste("Processing ", st, "\n", sep = ""))
    sub.data <- data[data$SS == st, ]
    if (metric != "SorCorr") {
      sub.data <- sub.data[sub.data$Measurement > 0, ]
    }
    sites.matrix <- matrix(nrow = length(unique(sub.data$SSBS)), 
                           ncol = length(unique(sub.data$SSBS)))
    i1 <- 1
    for (s1 in unique(sub.data$SSBS)) {
      i2 <- 1
      for (s2 in unique(sub.data$SSBS)) {
        if (metric == "Sor") {
          u <- length(union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                          s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                             s2]))
          i <- length(intersect(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                              s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                 s2]))
          sor <- (2 * i)/((2 * i) + (u - i))
        }
        else if (metric == "SorAbd") {
          u <- sum(sub.data$Measurement[(sub.data$SSBS == 
                                           s1) & (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                       s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                          s2]))])
          v <- sum(sub.data$Measurement[(sub.data$SSBS == 
                                           s2) & (sub.data$Taxon_name_entered %in% union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                       s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                          s2]))])
          sor <- (2 * u * v)/(u + v)
        }
        else if (metric == "SorCorr") {
          n <- sum(sub.data$Measurement[sub.data$SSBS == 
                                          s1])
          m <- sum(sub.data$Measurement[sub.data$SSBS == 
                                          s2])
          if ((n > 0) & (m > 0)) {
            xi <- sub.data$Measurement[sub.data$SSBS == 
                                         s1][(match(union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                        s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                           s2]), sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                               s1]))]
            yi <- sub.data$Measurement[sub.data$SSBS == 
                                         s2][(match(union(sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                        s1], sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                           s2]), sub.data$Taxon_name_entered[sub.data$SSBS == 
                                                                                                                                                               s2]))]
            xi[is.na(xi)] <- 0
            yi[is.na(yi)] <- 0
            f1. <- length(which((xi == 1) & (yi > 0)))
            f2. <- max(1, length(which((xi == 2) & (yi > 
                                                      0))))
            f.1 <- length(which((xi > 0) & (yi == 1)))
            f.2 <- max(1, length(which((xi > 0) & (yi == 
                                                     2))))
            p1 <- sum(xi[yi > 0]/n)
            p2 <- ((m - 1)/m) * (f.1/(2 * f.2))
            p3 <- sum(xi[yi == 1]/n)
            u <- min(1, p1 + p2 * p3)
            q1 <- sum(yi[xi > 0]/m)
            q2 <- ((n - 1)/n) * (f1./(2 * f2.))
            q3 <- sum(yi[xi == 1]/m)
            v <- min(1, q1 + q2 * q3)
            if ((u > 0) & (v > 0)) {
              sor <- (2 * u * v)/(u + v)
            }
            else {
              sor <- 0
            }
          }
          else {
            sor <- 0
          }
        }
        else {
          stop("Error: specfied dissimilarity metric is not supported")
        }
        if (s1 != s2) 
          sites.matrix[i1, i2] <- sor
        i2 <- i2 + 1
      }
      i1 <- i1 + 1
    }
    site.lu <- paste(sub.data$LandUse)[match(unique(sub.data$SSBS), 
                                             sub.data$SSBS)]
    lu.matrix <- matrix(nrow = length(unique(site.lu)), ncol = length(unique(site.lu)))
    for (lu1 in 1:length(unique(site.lu))) {
      for (lu2 in 1:length(unique(site.lu))) {
        lu.matrix[lu1, lu2] <- mean(sites.matrix[which(site.lu == 
                                                         unique(site.lu)[lu1]), which(site.lu == unique(site.lu)[lu2])], 
                                    na.rm = T)
      }
    }
    if (("Primary Vegetation" %in% site.lu) & !(TRUE %in% 
                                                  is.na(lu.matrix))) {
      used.studies <- used.studies + 1
      lu.matrix <- lu.matrix/lu.matrix[which(unique(site.lu) == 
                                               "Primary Vegetation"), which(unique(site.lu) == 
                                                                              "Primary Vegetation")]
      count.matrix[match(unique(site.lu), all.lu), match(unique(site.lu), 
                                                         all.lu)] <- count.matrix[match(unique(site.lu), 
                                                                                        all.lu), match(unique(site.lu), all.lu)] + 1
      sum.matrix[match(unique(site.lu), all.lu), match(unique(site.lu), 
                                                       all.lu)] <- sum.matrix[match(unique(site.lu), 
                                                                                    all.lu), match(unique(site.lu), all.lu)] + lu.matrix
      if (TRUE %in% is.na(sum.matrix)) 
        stop("NAs in dissimilarity matrix")
    }
  }
  all.results <- sum.matrix/count.matrix
  all.results <- data.frame(all.results)
  names(all.results) <- all.lu
  row.names(all.results) <- all.lu
  all.results[is.na(all.results)] <- NA
  cat(paste("Considered ", used.studies, " of ", all.studies.count, 
            "\n", sep = ""))
  return(list(cd = all.results, n = count.matrix))
}
# Also insert grouping for predictsdata
d <- data.frame(from=unique(predictsdata$Study_common_taxon),
                to=c("Invertebrates","Invertebrates","Invertebrates","Birds","Reptiles","Mammals","Other",
                     "Invertebrates","Plants","Amphibia","Plants","Plants","Plants","Plants","Mammals","Other","Plants","Amphibia",
                     "Mammals","Invertebrates")) # for matching
predictsdata$Grouping <- NA
predictsdata$Grouping <-  d$to[match(predictsdata$Study_common_taxon,d$from)]
rm(d)

# average pairwise comparison between primary and cropland
pd2 <- subset(predictsdata,Grouping=="Birds")# Only Birds
#pd3 <- subset(predictsdata,Grouping=="Plants")# Only Plants

pd <- predictsdata
compdis <- modCompDiss(pd,metric = "Sor",T)
compdis_b <- modCompDiss(pd2,metric = "Sor",T)
#compdis_p <- modCompDiss(pd3,metric = "Sor",T)

sTmerge <- subset(sites,select=c("SiteName","Transect"));sTmerge$SS <- ifelse(sTmerge$Transect=="Taita","MJ-Taita","MJ-Kilimanjaro");sTmerge$Transect <- NULL;names(sTmerge) <- c("SSBS","SS")
# Have to reshape fielddata to match PREDICTS data layout
s <- melt(species);names(s) <- c("SSBS","Taxon_name_entered","Measurement") # wide to long
# Subset
s$Diversity_metric_type <- "Abundance"
s$Diversity_metric_unit <- "individuals"
s <- plyr::join(s,sTmerge,by="SSBS",match = "first")
s$Predominant_habitat <- sites$PREDICTS.LU[match(s$SSBS,sites$SiteName)]
#s <- subset(s,Predominant_habitat!="Urban") # To be consistent with the birds only data
compdis_f <- modCompDiss(s,metric = "Sor")

# The same for Dickens Study
#sTmerge <- subset(sites,select=c("SiteName","Transect"));sTmerge$SS <- ifelse(sTmerge$Transect=="Taita","DO-Taita","DO-Kilimanjaro");sTmerge$Transect <- NULL;names(sTmerge) <- c("SSBS","SS")
#s <- melt(as.matrix(plantdata));names(s) <- c("SSBS","Taxon_name_entered","Measurement") # wide to long
#s$Diversity_metric_type <- "Abundance"
#s$Diversity_metric_unit <- "individuals"
#s$SS <- "DO-Fieldwork"
##s <- plyr::join(s,sTmerge,by="SSBS",match = "first") # Might not do due to missing PF
#s$Predominant_habitat <- sites$PREDICTS.LU[match(s$SSBS,sites$SiteName)]
#compdis_f2 <- modCompDiss(s,metric = "Sor")

N = data.frame(id=rownames(compdis$cd),n=compdis$n[1,]);N <- N[order(ordered(N$id,ord)),]
N_B = c(max(compdis$n),max(compdis_b$n) )
#N_P = c(max(compdis$n),max(compdis_p$n) )

# Make plot
rect_colours<-c("#00AE00","#94BD5E", "#006B6B", "#E6E64C","#808080")
rect_lab <- c("PV","SV","PL","CL","UR")
d <- rbind( 
            data.frame(compdis$cd[1,],typ="PREDICTS",fill=paste0("Full data")),
            data.frame(cbind(compdis_b$cd[1,],Urban=NA),typ="PREDICTS",fill=paste0("Birds only")),
#            data.frame(compdis_p$cd[1,],typ="PREDICTS",fill=paste0("Plants only (N=",N_P[2],")")),
            # Fieldwork            
            data.frame(rev(compdis_f$cd[nrow(compdis_f$cd),]),typ="Fieldwork",fill=paste0("Birds only")),
            data.frame(rev(compdis_f$cd[nrow(compdis_f$cd),]),typ="Fieldwork",fill=paste0("Full data"))

#           ,data.frame(cbind(rev(compdis_f2$cd[nrow(compdis_f2$cd),]),Urban=NA),typ="Fieldwork",fill=paste0("Plants only (N=",N_P[2],")"))
)
            
d <- melt(d,id.vars = c("typ","fill"))
#d <- d[which(!d$variable=="Urban"),] # Kickout urban for now
d$variable <- factor(d$variable,c("Primary.Vegetation","Secondary.Vegetation","Plantation.forest","Cropland","Urban"))
#d$col<- rep(rect_colours,each = 2)
# Rename Levels of typ for facet headers
levels(d$typ) <- c("Broad-scale data","Fieldwork")
rect_colours<-c("#00AE00","#94BD5E", "#006B6B", "#E6E64C","#808080")

ggthemr("fresh",layout = "clean",type = "outer",spacing = 0.5)
#ggthemr("greyscale",layout = "clean")
#theme_set(theme_classic(base_size=16,base_family = "sans"))
g <- ggplot(d,aes(x=variable,y=value,group=fill,fill=variable))
g <- g + geom_hline(yintercept=c(.5,1), lty=2,size=.5,alpha=.6)
g <- g + geom_bar(stat="identity",position="dodge",color="black",width=.75,alpha=.95) #+ scale_fill_manual(values = rect_colours)
g <- g + facet_grid(fill~typ)
g <- g + scale_x_discrete(labels=rect_lab)
g <- g + scale_y_continuous(labels = percent_format(),expand=c(0,0),limits=c(0,1.1))
g <- g + xlab("") + ylab("Assemblage similarity with baseline") # No labels  
g <- g + scale_fill_manual(values=rect_colours,guide = guide_legend(direction = "horizontal",title.hjust = 0.5, title.position = "top",title="Land-use",
                                              label.position="bottom", label.hjust = 0.5, label.vjust = 0.5,
                                              label.theme = element_text(angle = 90,colour = "white"))) 
g <- g + theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(),legend.position="none") # Remove legend and x-axis stuff
g <- g + ggtitle("Similarity between assemblages in both datasets")
g
#ggsave("../Conferences/Conservation Symposium/CommOverlapBirdsOnly.png",plot=g,scale=1.1,dpi=400)
ggsave("OutputComparison/Comparison_CommunityOverlap_pub.png",plot=g,dpi=400)
#FOR THESIS
ggsave("../000Thesis_writeup/WriteUpLatex/gfx/Figure6.pdf",plot=g+theme(legend.position="none"),scale=1.1,dpi=400)

# Transparent
t <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "transparent",colour = NA),
           axis.text = element_text(colour="white"),axis.text.x =element_text(colour="white"),axis.text.y =element_text(colour="white"),axis.ticks = element_line(colour = "white"), axis.line = element_line(colour="white"),axis.title = element_text(colour="white"),
           legend.key = element_rect(fill=NA,colour = NA),legend.background = element_rect(fill = NA,colour=NA),legend.text = element_text(color="white"),
           title = element_text(colour="white"),
           strip.background = element_rect(fill=NA,color=NA), strip.text = element_text(color="white"),
           panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
           text = element_text(color="white")
) 
ggsave("../001Thesis_Defense/assets/img/Figure6_trans.png",plot=g+t,scale = 1.1,units = "mm",dpi=400,bg="transparent")



#### Do both data have similarities in remote-sensing data- boxplot ####
sites$logpop <- log1p(sites$pop)
sites.pred$logpop <- log1p(sites.pred$pop)
sites$Type <- sites$Transect;sites.pred$Type <- "Africa-wide"
sites$Species_richness <- sites$spec
sites$SS <- sites$Transect
sites <- subset(sites,Grouping=="Aves") # Get red of Dickens data
# Do Figure for all Land-Use Types. Grouped by Fieldwork | PREDICTS
rm(full,full2)
preds <- c("SS","Species_richness","logabund","PREDICTS.LU","PREDICTS.LUI","logpop","meanNDVI","yieldNDVI","FC2000","Type","Grouping","yield.ndvi.corr","HYDE_Change")
full <- rbind(sites.pred[,preds],sites[,preds])
full$Type <- ordered(full$Type,rev(c("Africa-wide","Taita","Kilimanjaro")))
full$Type2 <- as.factor(ifelse(full$Type=="Africa-wide","Africa-wide","Fieldwork"))
full <- full[which(!is.na(full$PREDICTS.LU)),] # Kickout NA
full2 <- melt(full,na.rm = T,id.vars = c("Type","PREDICTS.LU","PREDICTS.LUI"),measure.vars = c("logpop","yield.ndvi.corr","FC2000","meanNDVI") )

levels(full2$variable) <- c("(log) Population\n Density","yield NDVI corrected","Forest Cover\n (in 2000)","mean NDVI")
label <- c("PV","SV","PL","CL","UR")

#ggthemr("earth",layout = "clean",type = "outer",spacing = 0)
#ggthemr("fresh",layout = "clean",type = "outer",spacing = 0.5)
#ggthemr("greyscale",layout = "clean")
theme_set(theme_classic(base_size=12,base_family = "sans"))
p <- ggplot(full2,aes(x=PREDICTS.LU,y=value,fill=Type))
p <- p + geom_boxplot(outlier.colour=NA,alpha=.5)
p <- p + facet_wrap(~variable,nrow = 2,drop = T,scales = "free_y")
p <- p + scale_x_discrete(labels=label)
#p <- p + theme(strip.text.x = element_text(size=10),strip.background = element_rect(colour="black", fill="white"),
#               axis.ticks.x=element_blank())
p <- p + theme(legend.position="bottom")
p <- p + scale_fill_brewer(palette = "Dark2",guide=guide_legend(direction = "horizontal",title.position = "left",title="Sites:",
                                                                            label.position="bottom", label.hjust = 0.5, label.vjust = 0.5))
p <- p + theme(plot.margin=unit(c(1, 1, -0.5, 0.5), units="line"))
p <- p + ylab("Value range of measured predictor") + xlab("")
p <- p + theme(panel.margin=unit(0,"lines"))
p 
ggsave("OutputComparison/AllPredictors.png",plot=p + scale_fill_manual(values=c("Grey40","Grey90","Black")) ,dpi=400)
#ggsave("../Conferences/Conservation Symposium/AllPredictors.png",scale = 1.1,plot=p,dpi=400)
#FOR THESIS
ggsave("../000Thesis_writeup/WriteUpLatex/gfx/Figure3.pdf",plot=p,scale=1.1,dpi=400)

# Transparent
t <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "transparent",colour = NA),
           axis.text = element_text(colour="white"),axis.text.x =element_text(colour="white"),axis.text.y =element_text(colour="white"),axis.ticks = element_line(colour = "white"), axis.line = element_line(colour="white"),axis.title = element_text(colour="white"),
           legend.key = element_rect(fill=NA,colour = NA),legend.background = element_rect(fill = NA,colour=NA),legend.text = element_text(colour="white"),
           title = element_text(colour="white"),
           strip.background = element_rect(fill=NA,color=NA), strip.text = element_text(color="white"),
           panel.grid.major = element_blank(),panel.grid.minor = element_blank()
) 
ggsave("../001Thesis_Defense/assets/img/Figure3_trans.png",plot=p+t,scale = 1.1,units = "mm",dpi=400,bg="transparent")


ddply(sites,.(Transect),summarise,Avg=mean(HYDE_Change,na.rm=T),SE=sem(HYDE_Change))
ddply(sites.pred,.(Source_ID),summarise,Avg=mean(HYDE_Change,na.rm=T),SE=sem(HYDE_Change))


summarySE(sites,"HYDE_Change","Transect",na.rm=T)
summarySE(sites.pred,"HYDE_Change",na.rm=T)

# Now build an Intensity score and plot weighted average for LU
d <- data.frame(PREDICTS.LUI=unique(full2$PREDICTS.LUI),sc=c(1,2,3,0))
full2 <- join(full2,d)

wm <- ddply(full2,.(PREDICTS.LU,Type2,variable),summarise,
            wm = weighted.mean(value,sc))

library(ggthemr)
ggthemr("greyscale",layout = "clean")
#theme_set(theme_classic(base_size=12,base_family = "sans"))
p <- ggplot(wm,aes(x=PREDICTS.LU,y=wm,color=Type2))
p <- p + geom_point()
p <- p + facet_wrap(~variable,nrow = 2,drop = T,scales = "free_y")
p <- p + scale_x_discrete(labels=label)
p <- p + ylab("Value range of measured predictor") + xlab("")
p <- p + theme(panel.margin=unit(0,"lines"))
p


# Plot all sites abundace / SR against yieldEVI and meanEVI and logpop
# smoothed line plot
plot(Species_richness~yieldEVI,data=subset(full,Grouping=="Birds"),col=Type)
plot(logabund~yieldEVI,data=subset(full,Grouping=="Birds"),col=Type)
plot(logabund~meanEVI,data=subset(full,Grouping=="Birds"),col=Type)
plot(Species_richness~meanEVI,data=subset(full,Grouping=="Birds"),col=Type)
plot(Species_richness~logpop,data=subset(full,Grouping=="Birds"),col=Type)


plot((full$Species_richness)~full$meanEVI,col=full$Type)
plot((full$logabund)~full$meanEVI,col=full$Type)
plot((full$Species_richness)~full$logpop,col=full$Type)

ggplot(full) + geom_density(aes(x = new.yield.evi,fill = Type),alpha=.3)
ggplot(full) + geom_density(aes(x = logpop,fill = Type),alpha=.3)

g <- ggplot(full,aes(x = new.yield.evi,y=logpop))
g <- g + geom_density2d(aes(col = Type),size=1.2)
g <- g + scale_color_brewer(palette = "Dark2") + scale_size_area()
g <- g + xlim(c(0, .2)) + ylim(c(-4.5, 3))
g

pairs(full[,c("new.yield.evi","new.mean.evi","logpop")])


#### Custom Table of Authors ####
# Just make a custom table taking the author/study list from Lawrence
studies <- read.csv("Data/martin-sources.csv",header=T)
studies$Insightly_restrictions <- NULL; studies$Insightly_category <- NULL # Kickout insightly stuff
studies <- studies[which(studies$Source_ID%in%sites.pred$Source_ID),] # Get only those studies that I analysed

a = ddply(sites.pred,("Source_ID"),summarise,
        NrSites = length(SSS),
        LandUseClasses = paste(as.character(unique(PREDICTS.LU)),collapse =  ","),
        Taxon = paste(as.character(unique(Grouping)),collapse = ",")        
      )
b <- ddply(predictsdata,("Source_ID"),summarise,
           SpeciesRichnes = length(unique(Best_guess_binomial))
           )
aa <- join(a,b)
studies <- join(studies,aa);rm(a,b,aa)

studies$Source_ID <- NULL
#studies$First_author_surname <- NULL
studies$Publication_type <- NULL
studies <- studies[order(studies$Year_of_publication),]

write.csv(studies,"../000Thesis_writeup/Chapter1/SI_ContributedStudies.csv",row.names=F)


#### Other Appendix figures ####
source("../../Work/Analysis_Fruitgrowth/ggplot2_imp.R")
ds <- sites.pred[which(!is.na(sites.pred$LandUse)),]
dfc <- marfunky::summarySEwithin(ds, measurevar="logabund",na.rm=T, withinvars=c("UN_subregion","LandUse"))
dfc$LandUse <- factor(dfc$LandUse,levels = ord)
theme_set(theme_classic(base_size=12,base_family = "sans"))
g1 <- ggplot(dfc,aes(x=LandUse,y=logabund,fill=UN_subregion,group=UN_subregion))
g1 <- g1 + geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=logabund-se, ymax=logabund+se,width=0.2))
g1 <- g1 + geom_bar(stat="identity", position="dodge") # WORKAROUND: once with no colour, for the legend...
g1 <- g1 + geom_bar(stat="identity", position="dodge", colour="black",show_guide=FALSE)#
g1 <- g1 + theme(axis.ticks.x = element_blank()) # Remove x-Axis ticks
g1 <- g1 + xlab("") + ylab(expression(log("Abundance"))) # No labels  
g1 <- g1 + ggtitle("Geographic representation \n of broad-scale abundance data")
g1 <- g1 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g1 <- g1 + scale_fill_brewer(palette="Spectral",guide = guide_legend(title.position = "top",title="UN Subregion",
                                                                   label.position="right"))

ggsave(filename="../000Thesis_writeup/WriteUpLatex/gfx/FigureAppendixAdd1.pdf",plot=g1,scale=1.1,dpi=400)

# Transparent
t <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "transparent",colour = NA),
           axis.text = element_text(colour="white"),axis.text.x =element_text(colour="white"),axis.text.y =element_text(colour="white"),axis.ticks = element_line(colour = "white"), axis.line = element_line(colour="white"),axis.title = element_text(colour="white"),
           legend.key = element_rect(fill=NA,colour = NA),legend.background = element_rect(fill = NA,colour=NA),legend.text = element_text(colour="white"),
           panel.grid.minor = element_blank(), title = element_text(colour="white")
)
ggsave(filename="../001Thesis_Defense/assets/img/SI1_trans.png",plot=g1+t,scale=1.1,dpi=400,bg = "transparent")

library(dplyr)
x <- ds %>% group_by(PREDICTS.LU,Grouping) %>%
  summarise(n = n()) %>% group_by(PREDICTS.LU) %>%  mutate(prop = n /sum(n)) %>%
  group_by(PREDICTS.LU,Grouping) %>% arrange(desc(prop)) %>%
  as.data.frame

g2 <- ggplot(x ,aes(x=PREDICTS.LU,y=prop,fill=Grouping))
g2 <- g2 + geom_bar(stat="identity")
g2 <- g2 + xlab("") + ylab("Proportion of Sites")
g2 <- g2 + scale_y_continuous(labels=percent_format())
g2 <- g2 + theme(axis.ticks.x = element_blank())
g2 <- g2 + ggtitle("Taxonomic representation \n in African PREDICTS Data")
g2 <- g2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g2 <- g2 + scale_fill_brewer(palette="Dark2",guide = guide_legend(title.position = "top",title="Taxonomic affilation",
                                                                   label.position="right"))
g2
ggsave(filename="../000Thesis_writeup/WriteUpLatex/gfx/FigureAppendixAdd2.pdf",plot=g2,scale=1.1,dpi=400)

# Transparent
ggsave(filename="../001Thesis_Defense/assets/img/SI2_trans.png",plot=g2+t,scale=1.1,dpi=400,bg = "transparent")
