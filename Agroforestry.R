##### Package Loading----
library(marfunky)
standardPackages()
library(roquefort)
library(yarg)
library(rgdal)
library(raster)
# Load PREDICTS and Sites default
sem <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

# Load Agroforestry
ag <- read.csv("Agroforestry.csv",header=T,sep=";")
sub.sites.pred <- merge(sub.sites.pred,ag,by=c("SS","Site_name","LandUse"),all.x = T)

#### Sample From Global fieldSizes -----
crop_f <- raster("../GIS/IASA_CROPLAND/CroplandSize_Africa.tif")
cropp <- raster("../GIS/IASA_CROPLAND/CroplandProportion.tif")

sub.sites.pred$FieldSize <- with(sub.sites.pred,{
  FieldSize <- extract(crop_f, SpatialPoints(cbind(Longitude.x,Latitude.x) ))
        })
sub.sites.pred$CropPop <- with(sub.sites.pred,{
  FieldSize <- extract(cropp, SpatialPoints(cbind(Longitude.x,Latitude.x) ))
})
# Field
sites$FieldSize <- with(sites,{
  FieldSize <- extract(crop_f, SpatialPoints(cbind(Long,Lat) ))
})
sites$CropPop <- with(sites,{
  FieldSize <- extract(cropp, SpatialPoints(cbind(Long,Lat) ))
})

#### Compare Predictors of Homegardens ####
sites$logpop <- log1p(sites$pop)
sub.sites.pred$logpop <- log1p(sub.sites.pred$pop)
sites$Type <- sites$Transect;sub.sites.pred$Type <- "Africa-wide"
sites$Species_richness <- sites$spec
sites$SS <- sites$Transect
sites$Agroforestry <- sites$Homegarden == 1
sites$AF_Type <- ifelse(sites$Homegarden==1,"Homegarden","")
sites <- subset(sites,Grouping=="Aves") # Get red of Dickens data
# Do Figure for all Land-Use Types. Grouped by Fieldwork | PREDICTS
rm(full,full2)

preds <- c("SS","Species_richness","logabund","PREDICTS.LU","PREDICTS.LUI","logpop","meanNDVI","yieldNDVI","FC2000","Type","Grouping","yield.ndvi.corr","HYDE_Change",
           "Agroforestry","AF_Type","FieldSize","CropPop")
full <- rbind(sub.sites.pred[,preds],sites[,preds])
full$Type <- ordered(full$Type,rev(c("Africa-wide","Taita","Kilimanjaro")))
full$Type2 <- as.factor(ifelse(full$Type=="Africa-wide","Africa-wide","Fieldwork"))
full <- full[which(!is.na(full$PREDICTS.LU)),] # Kickout NA
full2 <- melt(full,na.rm = T,id.vars = c("Type2","PREDICTS.LU","PREDICTS.LUI"),measure.vars = c("yield.ndvi.corr","meanNDVI","FieldSize","CropPop") )

levels(full2$variable) <- c("yield NDVI corrected","mean NDVI","FieldSize","Cropland Proportion")
label <- c("PV","SV","PL","CL","UR")

library(ggthemr)
ggthemr("fresh",layout = "clean",type = "outer",spacing = 0.5)
#ggthemr("greyscale",layout = "clean")
p <- ggplot(full2,aes(x=Type2,y=value,fill=PREDICTS.LU))
p <- p + geom_boxplot(outlier.colour=NA,alpha=.5)
p <- p + facet_wrap(~variable,nrow = 2,drop = T,scales = "free_y")
#p <- p + theme(strip.text.x = element_text(size=10),strip.background = element_rect(colour="black", fill="white"),
#               axis.ticks.x=element_blank())
p <- p + theme(legend.position="bottom")
p <- p + scale_fill_brewer(palette = "Dark2",guide=guide_legend(direction = "horizontal",title.position = "left",title="Sites:",
                                                                label.position="bottom", label.hjust = 0.5, label.vjust = 0.5))
p <- p + theme(plot.margin=unit(c(1, 1, -0.5, 0.5), units="line"))
p <- p + ylab("Value range of measured predictor") + xlab("")
p <- p + theme(panel.margin=unit(0,"lines"))
p
#ggsave("../000Thesis_writeup/WriteUpLatex/gfx/Figure3.pdf",plot=p,scale=1.1,dpi=400)


# Show Agroforestry only
f <- subset(full,(PREDICTS.LU=="Cropland"&PREDICTS.LUI=="Minimal use")|(Agroforestry==T&AF_Type=="Homegarden"))
#f <- subset(full,Agroforestry==T&AF_Type=="Homegarden")
f$PREDICTS.LU <- droplevels(f$PREDICTS.LU)
# Make Index
f$AF_Type <- as.character(f$AF_Type)
f$AF_Type[f$Agroforestry==F] <- paste0(f$PREDICTS.LU[f$Agroforestry==F],"\n(",f$PREDICTS.LUI[f$Agroforestry==F],")")
f$AF_Type[f$AF_Type=="Homegarden"] <- "Tropical Homegarden"
f$AF_Type <- as.factor(f$AF_Type)
full2 <- melt(f,na.rm = T,id.vars = c("AF_Type","Type2"),measure.vars = c("yield.ndvi.corr","meanNDVI","FieldSize","CropPop") )
levels(full2$variable) <- c("yield NDVI corrected","mean NDVI","Field Size (10-40 ordinal) ","Cropland proportion (%)")

#ddply(f,.(Type2,AF_Type),summarise,N=n())

library(ggthemr)
ggthemr("fresh",layout = "clean",type = "outer",spacing = 0.5)
#ggthemr("greyscale",layout = "clean")
p <- ggplot(full2,aes(x=AF_Type,y=value,fill=Type2))
p <- p + geom_boxplot(outlier.colour=NA,alpha=.5)
p <- p + facet_wrap(~variable,nrow = 2,drop = T,scales = "free_y")
#p <- p + theme(strip.text.x = element_text(size=10),strip.background = element_rect(colour="black", fill="white"),
#               axis.ticks.x=element_blank())
p <- p + theme(legend.position="bottom")
p <- p + scale_fill_brewer(palette = "Set1",guide=guide_legend(direction = "horizontal",title.position = "left",title="Sites:",
                                                                label.position="bottom", label.hjust = 0.5, label.vjust = 0.5))
p <- p + theme(plot.margin=unit(c(1, 1, -0.5, 0.5), units="line"))
p <- p + ylab("Value range of measured predictor") + xlab("") + ggtitle("Delineation of Agroforestry (Tropical Homegardens)")
p <- p + theme(panel.margin=unit(0,"lines"))
p
ggsave("../000Thesis_writeup/WriteUpLatex/gfx/Figure7.pdf",plot=p,scale=1.1,dpi=400)


#### MODEL Comparison ####
# Integrate Agroforestry
ord <- c("Primary Vegetation", "Secondary Vegetation","Plantation forest","Pasture","Agroforestry", "Cropland","Urban")

sub.sites.pred$PREDICTS.LU <- as.character(sub.sites.pred$PREDICTS.LU)
sub.sites.pred$PREDICTS.LU[with(sub.sites.pred,{Agroforestry==T&AF_Type=="Homegarden"})] <- "Agroforestry"
sub.sites.pred$PREDICTS.LU <- factor(sub.sites.pred$PREDICTS.LU,levels = ord)
# For Field sites
sites$PREDICTS.LU <- as.character(sites$PREDICTS.LU)
sites$PREDICTS.LU[with(sites,{Agroforestry==T&AF_Type=="Homegarden"})] <- "Agroforestry"
sites$PREDICTS.LU <- factor(sites$PREDICTS.LU,levels = ord)

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

rc <- CompareRandoms(sub.sites.pred,"logabund","gaussian",fixedFactors=c("PREDICTS.LU"),siteRandom=F)
pfit_lu_simple <- glmer(as.formula(paste("logabund ~  PREDICTS.LU +",rc$best.random)),data=sub.sites.pred,family=gaussian,REML=F)
rc <- CompareRandoms(sub.sites.pred,"Species_richness","poisson",fixedFactors="PREDICTS.LU",siteRandom=T)
pfit_lu_simple_rich <- glmer(as.formula(paste("Species_richness ~ PREDICTS.LU + ",rc$best.random)),data=sub.sites.pred,family=poisson(link=log),REML=F)
dfc <- makeRatios(pfit_lu_simple,"PREDICTS.LU",T,F);dfc$data <- "All";dfc$type <- "log(Abundance)";dfc$Grouping <- NA
dfc_rich <- makeRatios(pfit_lu_simple_rich,"PREDICTS.LU",T,F);dfc_rich$data <- "All";dfc_rich$type <- "Species Richness";dfc_rich$Grouping <- NA

saveRDS(sub.sites.pred,"~/sub.sites.pred.rds")

# Site Data
# Build indices for each Grouping
sA <- subset(sites,Grouping=="Aves")

detach("package:dplyr", unload=TRUE)
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

sA$val <- sA$logabund
ofA <- sitemakeRatios(sA,"PREDICTS.LU",exp=T,"Birds")
ofA$type <- "log(Abundance)"

sA$val <- sA$spec
ofA_r <- sitemakeRatios(sA,"PREDICTS.LU",exp=F,"Birds")
ofA_r$type <- "Species Richness"

full <- rbind(ofA,ofA_r,dfc,dfc_rich)#,ofP,ofP_r)
full$LandUse <- ordered(full$LandUse,ord)
full$Transect <- ordered(full$Transect,levels=c("Kilimanjaro","PREDICTS","Taita"))

# Pub plot
rect_lab <- c("PV","SV","PL","AF","CL","UR")
xmin = c(0.5,1.5,2.5,3.5,4.5,5.5);xmax=c(1.5,2.5,3.5,4.5,5.5,6.5)
library(ggthemr)
#ggthemr("greyscale",layout = "clean")
ggthemr("fresh",layout="clean",type="outer",spacing=0)
g <- ggplot(full,aes(x=LandUse,y=eff,shape=Transect,group=Transect))
g <- g + geom_hline(yintercept=100, lty=2,size=.5)
g <- g + geom_vline(xintercept=xmax[-6],size=.2,linetype="dotted") + geom_vline(xintercept=0)# and vline seperation
g <- g + geom_pointrange(aes(ymax=se.high,ymin=se.low),size=.9,position=position_dodge(width=1) )
g <- g + facet_wrap(~type,scales = "free_y",as.table = F,drop = T) #+ theme(strip.text.x = element_text(size=16), strip.background = element_rect(colour="black"))
g <- g + scale_y_continuous(expand=c(0,1),breaks=pretty_breaks()) + coord_cartesian()
g <- g + scale_x_discrete(labels=rect_lab)
g <- g + labs(y="Percentage of Baseline",x="")
g <- g + scale_shape_manual(values=c(15,16,17),labels = c("Kilimanjaro","Broad-scale","Taita"),
                            guide = guide_legend(title = "Dataset:",direction = "horizontal", title.position = "left",
                                                 label.position="right", label.hjust = 0.5, label.vjust = 0.5,nrow=1)) 
g <- g + theme(legend.position="bottom",axis.title.x=element_blank(),axis.ticks.x=element_blank())
g
ggsave("../000Thesis_writeup/WriteUpLatex/gfx/Figure8_Agroforestry.pdf",plot=g,scale = 1.1,units = "mm",dpi=400)
