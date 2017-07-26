#### Preperation and standard Package Loading ####
# Packages
library(raster)
library(gdata)
library(BiodiversityR)
library(vegan)
library(rgdal)
library(ggplot2)
library(scales)
library(reshape2)
rm(list=ls())
par.ori <- par(no.readonly=T)
source("~/Dokumente/Studium/Kopenhagen/Masters/FieldData/PrepareFiles.R")
ori.dir = "~/Dokumente/Studium/Kopenhagen/Masters/Analysis_PREDICTS/"
setwd(ori.dir)
if(length(list.dirs("OutputDiversityAnalyis"))==0) dir.create("OutputDiversityAnalyis")
output.dir = "OutputDiversityAnalyis/"
# Get data

#### Basic diagnostic plots ####
# Most common species and some plots
head(sort(colSums(species),decreasing=T))
plot(sites$abund~sites$Elev,col=sites$Zone)

plot(sites$abund~sites$Elev,col=sites$Transect)
plot(sites$SpecCV~sites$Elev,col=sites$Transect)

summary(lm(log(sites$abund)~sites$Elev))
summary(lm(log(sites$shan)~sites$Elev))

barplot(table(sites$Zone))
boxplot(sites$Elev~sites$PREDICTS.LU,las=2)
boxplot(sites$Elev~sites$Transect,las=2)


# Elevational distribution of Landuse classes in both transects
theme_set(theme_classic(base_size=12,base_family = "sans"))
g <- ggplot(sites, aes(x=PREDICTS.LU, y=Elev,color=Transect, fill=Transect))
g <- g + geom_boxplot(outlier.size=0) + scale_color_grey(start=.6,end=.7) # WORKAROUND: once with no colour, for the legend...
g <- g + scale_y_continuous("Altitude",expand=c(0,0),limit=c(500,max(sites$Elev)))
g <- g + theme(axis.ticks.x = element_blank()) # Remove x-Axis ticks
g <- g + xlab("") # No labels  
g <- g + ggtitle("Distribution of Land-Use classes accross Altitude")
g <- g + scale_fill_grey(start=.4,end=.5,guide = guide_legend(title.position = "top",label.position="right"))
g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
g
ggsave(filename=paste0(output.dir,"Landuse_ElevationDistribution.png"),plot=g,dpi=400)

# Abundance across LU and LUI
g <- ggplot(sites, aes(x=PREDICTS.LU, y=abund,color=Transect, fill=Transect))
g <- g + geom_boxplot(outlier.size=3,outlier.colour="grey60") + scale_color_grey(start=.6,end=.7) + coord_trans(y = "log10") + coord_flip()
g <- g + scale_fill_grey(start=.4,end=.5,guide = guide_legend(title.position = "top",label.position="right"))
g <- g + xlab("") + ylab("Log10(Abundance)") + ggtitle("Distribution of Land-Use classes accross Altitude")
g
ggsave(filename=paste0(output.dir,"GG_Landuse_Abundance.png"),plot=g,dpi=400)

dfc <- melt(sites, id=c("Transect", "PREDICTS.LU"),measure.vars="abund")

pdf(paste0(output.dir,"N_Landuse_Abundance.pdf"))
par(mfrow=c(2,1),pty="s",las=1, bty="l",pch=19, cex=1.0,cex.lab=0.5,cex.axis=0.7, tcl=-0.2, mar=c(2,3,1,0)+0.1, mgp=c(2.5,0.5,0))
boxplot(sites$abund~sites$PREDICTS.LU,las=1,main="Abundance",horizontal=T)
labels=aggregate(sites$abund,by=list(sites$PREDICTS.LU),FUN=length)
axis(4,labels$Group.1,paste("N=",labels$x),tick=F)
boxplot(sites$spec~sites$PREDICTS.LU,las=1,main="Species",horizontal=T)
axis(4,labels$Group.1,paste("N=",labels$x),tick=F)
par(par.ori)
dev.off()

# General collector accumulation per Zone
pdf(paste0(output.dir,"Altitude_CollectorSpecAccum.pdf"))
par(pty="s",las=1, bty="l",pch=19, cex=1.0,cex.lab=1,cex.axis=0.7, tcl=-0.2, mar=c(5,1,3,0)+0.1, mgp=c(2.5,0.5,0))
accumcomp(species,y=sites,factor="Zone",type="l",method="collector",label=F,legend=F,conditioned=T,permutations=100)
title("Sampling Altitude Zones")
legend("bottomright",legend=c("High","Mid","Low"),lty=1,col=c("Red","Blue","Green"),bty="n")
par(par.ori)
dev.off()

# Exact Accumulation per LandUse and Forest
pdf(paste0(output.dir,"SpecAccumLU.pdf"))
par(pty="s",las=1,xpd=T, bty="l",pch=19, cex=1.0,cex.lab=1.2,cex.axis=1, tcl=-0.2, mar=c(5,1,3,0)+0.1, mgp=c(2.5,0.5,0))
accumcomp(species,sites,factor="PREDICTS.LU",xlab="Individuals",labelit=F,type="l",method="rare",legend=F,conditioned=T,permutations=1000)
legend("bottomright",cex=1.5,legend=c("Cropland","Plantation","Primary forest","Secondary Vegetation","Urban"),lty=1,col=rainbow(5),bty="n")
title("Rarefaction curves")
par(par.ori)
dev.off()

par(pty="s",las=1, bty="l",pch=19, cex=1.0,cex.lab=1,cex.axis=0.7, tcl=-0.2, mar=c(5,1,3,0)+0.1, mgp=c(2.5,0.5,0))
accumcomp(species,sites,factor="Type",method="exact",rainbow=F,label=F,legend=F,conditioned=T,permutations=100)
title("Sampling Forest")
legend("bottomright",legend=c("Forest","Non-Forest"),lty=1,pch=c(1,2),bty="n")
par(par.ori)

#### Ordination ####
sk <- species[which(sites$Transect=="Kilimanjaro"),] # Make Local Subsets
sik <- subset(sites,Transect=="Kilimanjaro")
csik <- log1p(pointDistance(cbind(sik$Long,sik$Lat),lonlat=T,allpairs=T))
st <- species[which(sites$Transect=="Taita"),] 
sit <- subset(sites,Transect=="Taita")
csit <- log1p(pointDistance(cbind(sit$Long,sit$Lat),lonlat=T,allpairs=T))

# NMDS with log-transformed matrix using manhattan distances
d <- decostand(sk,method="log") # scale them to log
ord_k <- metaMDS(d,distance="manhattan",k=20,autotransform=F) # NMDS
d <- decostand(st,method="log") # scale them to log
ord_t <- metaMDS(d,distance="manhattan",k=20,autotransform=F) # NMDS

pdf(paste0(output.dir,"NMDS_Stressplots.pdf"))
par(mfrow=c(1,2),las=1, bty="l",pch=21, cex=1.0,cex.lab=1,cex.axis=0.7, tcl=-0.2, mgp=c(2,0.5,0))
stressplot(ord_k,main="Kilimanjaro",pch=".")
stressplot(ord_t,main="Taita",pch=".")
par(par.ori)
dev.off()


(fit <- envfit(ord_k~PREDICTS.LU,sik,permutations=999))
(fit2 <- envfit(ord~Landuse,sites,strata="Zone",permutations=999))
(fit3 <- envfit(ord~Zone,sites,permutations=999))

pdf(paste0(output.dir,"NMDS.pdf"))
fit <- envfit(ord_t~PREDICTS.LU,sit,permutations=999)
par(mfrow=c(1,2),pty="s",las=1, bty="l",pch=21, cex=1.0,cex.lab=1,cex.axis=0.7, tcl=-0.2, mar=c(3,3,3,0)+0.1, mgp=c(2,0.5,0))
plot(ord_t,type="p",display="sites",main="NMDS Taita")
#ordisurf(ord_t,sit$Elev,col=c("black"),lwd=2,nlevels=5,labcex=1,bubble=T,main="NMDS Taita")
ordihull(ord=ord_t,groups=sit$Zone,draw="polygon",col="lightblue",label=T,alpha=40,cex=.5)
plot(fit, labels=labels(fit),cex=0.7,p.max = 0.005,xpd=NA, col = "red") # Plot only significant Trends

fit <- envfit(ord_k~PREDICTS.LU,sik,permutations=999)
plot(ord_k,type="p",display="sites",main="NMDS Kilimanjaro")
#ordisurf(ord_t,sit$Elev,col=c("black"),lwd=2,nlevels=5,labcex=1,bubble=T,main="NMDS Taita")
ordihull(ord=ord_k,groups=sik$Zone,draw="polygon",col="lightblue",label=T,alpha=40,cex=.5)
plot(fit, labels=labels(fit),cex=0.7,p.max = 0.005,xpd=NA, col = "red") # Plot only significant Trends
par(par.ori)
dev.off()


# Specpool
(pool <- specpool(species,pool=sites$PREDICTS.LU))
par(mfrow=c(1,2),pty="s",las=1, bty="l",pch=19, cex=1.0,cex.lab=1,cex.axis=0.7, tcl=-0.2, mar=c(5,2,3,0)+0.1, mgp=c(2.5,0.5,0))
boxplot(specnumber(species) ~ PREDICTS.LU,data=sites,las=2)
boxplot(specnumber(species)/specpool2vect(pool) ~ PREDICTS.LU,data=sites,las=2)
par(par.ori)

p <- poolaccum(species,permutations=100)
plot(p,alpha=0.05,type="l")

x <- (nestedtemp(d,niter=20,names=T,weighted=T))
plot(x)
plot(x, kind="incid")
nestedchecker(d)
oecosimu(d,nestedchecker,"quasiswap") # significant non random nestedness pattern
nestedbetasor(species)

# Look at functcomp and dbFD in package FD

#### THESIS ####
# Make a plot of the 
hg <- sites[sites$Homegarden==1,]

qplot(hg$PREDICTS.LU,hg$Species_richness)
