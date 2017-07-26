#### Startup ####
#rm(list=ls())
par.ori <- par(no.readonly=T)
library(gdata)
library(vegan)
library(lme4)
library(MuMIn)
library(AICcmodavg)
library(influence.ME)
library(lubridate)
library(effects) 
library(lmerTest)
library(LMERConvenienceFunctions)
library(varComp)
library(spaMM)
library(coefplot2) # coefplot for lmerMod
#library(visreg) # For plotting interactions and regression surfaces
# Gam
library(mgcv)
library(gamm4)
library(ggplot2)
library(GGally)
library(rgdal)
library(raster)
library(gstat)
library(rgeos)
library(adehabitatMA)
library(reshape2)
library(picante)
# My own
library(marfunky)
standardPackages()
library(tidyr)
# Custom PREDICTS Packages
library(yarg)
library(roquefort)

ori.dir <- "/home/martin/Documents/Studium/Kopenhagen/Masters/Analysis_PREDICTS"
setwd(ori.dir)

# Get Field data covariates
source("../FieldData/PrepareFiles.R")
source("GIS_Loading.R")
sites <- field_data_load(sites)
sites$HYDE_Change <- AppendHYDE(sites$Long,sites$Lat,2014)

# Predicts data output from 28-07-2014
predictsdata <- readRDS("Data/martin-2014-07-28-16-47-26.rds")

predictsdata <- DropInvalidMetricsAndMethods(predictsdata)
predictsdata <- CorrectSamplingEffort(predictsdata)

# Calculate unique number of sites per land-use type
#pd <- subset(predictsdata,!(Predominant_habitat%in%c("Pasture")))
#pd <- predictsdata
#pd <- subset(pd,Source_ID!="GP1_2012__Strauch")
#pd <- subset(pd,year(sites.pred$Sample_start_earliest)>=2000 )
#ddply(pd,.(Predominant_habitat),summarise,N=length(unique(SSS)))

# Calculate number of species per Taxonomic group
#d <- data.frame(from=unique(pd$Study_common_taxon),to=c("Invertebrates","Invertebrates","Invertebrates","Birds","Reptiles","Mammals","Other","Invertebrates","Plants","Amphibia","Plants","Plants","Mammals","Other","Plants","Amphibia","Mammals","Invertebrates")) # for matching
#pd$Grouping <- NA
#pd$Grouping <-  d$to[match(pd$Study_common_taxon,d$from)]
#rm(d)
#ddply(pd,.(Grouping),summarise,N=length(unique(Species)))


# length(unique(pd$Source_ID))
#length(unique(pd$Species))
# length(unique(pd$SSS))
# rm(pd)
# length(sites$SiteName)
# ncol(species)

# diversity <- MergeSites(diversity)
sites.pred<-SiteMetrics(diversity=predictsdata,
                        extra.cols=c("SSB","SSBS","Ecoregion","Biome","Country","UN_subregion","Site_name",
                                     "Sampling_method","Study_common_taxon","Max_linear_extent",
                                     "Crop"),sites.are.unique=TRUE)

# Kick out unsuitable studies
sites.pred <- subset(sites.pred,Source_ID!="GP1_2012__Strauch")
sites.pred <- subset(sites.pred,!(Predominant_habitat%in%c("Pasture")))
sites.pred <- subset(sites.pred,year(sites.pred$Sample_start_earliest)>=2000)

# Reclassify LandUse values
sites.pred$LandUse<-paste(sites.pred$Predominant_habitat)
sites.pred$LandUse[which(sites.pred$LandUse=="Primary forest")]<-"Primary Vegetation"
sites.pred$LandUse[which(sites.pred$LandUse=="Primary non-forest")]<-"Primary Vegetation"
sites.pred$LandUse[which(sites.pred$LandUse=="Young secondary vegetation")]<-"Secondary Vegetation"
sites.pred$LandUse[which(sites.pred$LandUse=="Intermediate secondary vegetation")]<-"Secondary Vegetation"
sites.pred$LandUse[which(sites.pred$LandUse=="Mature secondary vegetation")]<-"Secondary Vegetation"
sites.pred$LandUse[which(sites.pred$LandUse=="Secondary vegetation (indeterminate age)")]<-"Secondary Vegetation"
sites.pred$LandUse[which(sites.pred$LandUse=="Secondary non-forest")]<-"Secondary Vegetation"
sites.pred$LandUse[which(sites.pred$LandUse=="Cannot decide")]<-NA
sites.pred$LandUse<-factor(sites.pred$LandUse)
sites.pred$LandUse<-relevel(sites.pred$LandUse,ref="Primary Vegetation")

sites.pred$LogSimpson<-log(sites.pred$Simpson_diversity)
sites.pred$logabund<-log1p(sites.pred$Total_abundance)


sites.pred$PREDICTS.LU <- sites.pred$LandUse
sites.pred$PREDICTS.LUI <- sites.pred$Use_intensity

# Split cropland minimal and light use
sites$CroplandSplit <- NA
sites$CroplandSplit[which(sites$PREDICTS.LU=="Cropland"&sites$PREDICTS.LUI=="Minimal use")] <- "Minimal use Cropland"
sites$CroplandSplit[which(sites$PREDICTS.LU=="Cropland"&sites$PREDICTS.LUI=="Light use")] <- "Light use Cropland"
sites$CroplandSplit[is.na(sites$CroplandSplit)] <- as.character(sites$PREDICTS.LU[is.na(sites$CroplandSplit)])

sites.pred$CroplandSplit <- NA
sites.pred$CroplandSplit[which(sites.pred$PREDICTS.LU=="Cropland"&sites.pred$PREDICTS.LUI=="Minimal use")] <- "Minimal use Cropland"
sites.pred$CroplandSplit[which(sites.pred$PREDICTS.LU=="Cropland"&sites.pred$PREDICTS.LUI=="Light use")] <- "Light use Cropland"
sites.pred$CroplandSplit[which(sites.pred$PREDICTS.LU=="Cropland"&sites.pred$PREDICTS.LUI=="Intense use")] <- "Intense use Cropland"
sites.pred$CroplandSplit[is.na(sites.pred$CroplandSplit)] <- as.character(sites.pred$PREDICTS.LU[is.na(sites.pred$CroplandSplit)])


# Re Leveling factor orders for regression
ord2 <- c("Primary Vegetation", "Secondary Vegetation","Plantation forest","Pasture", "Minimal use Cropland","Light use Cropland","Intense use Cropland","Urban")
ord <- c("Primary Vegetation", "Secondary Vegetation","Plantation forest","Pasture", "Cropland","Urban")
ordi <- c("Minimal use","Light use","Intense use","Cannot decide")

sites$PREDICTS.LU <- factor(sites$PREDICTS.LU, levels = ord)
sites$PREDICTS.LUI <- factor(sites$PREDICTS.LUI, levels = ordi)
sites.pred$PREDICTS.LU <- factor(sites.pred$PREDICTS.LU, levels = ord)
sites.pred$PREDICTS.LUI <- factor(sites.pred$PREDICTS.LUI, levels = ordi)

sites$CroplandSplit <- factor(sites$CroplandSplit, levels = ord2)
sites.pred$CroplandSplit <- factor(sites.pred$CroplandSplit, levels = ord2)

# Grouping
d <- data.frame(from=unique(sites.pred$Study_common_taxon),
                to=c("Invertebrates","Invertebrates","Invertebrates","Birds","Reptiles","Mammals","Other",
                     "Invertebrates","Plants","Amphibia","Plants","Plants","Mammals","Other","Plants","Amphibia",
                     "Mammals","Invertebrates")) # for matching
sites.pred$Grouping <- NA
sites.pred$Grouping <-  d$to[match(sites.pred$Study_common_taxon,d$from)]
rm(d)

#### Here####

# Load PREDICTS GIS data
print("Add PREDICTS predictors :)")
sites.pred <- predicts_data_load(sites.pred)
sites.pred$HYDE_Change <- AppendHYDE(sites.pred$Longitude,sites.pred$Latitude,year(sites.pred$Sample_start_earliest))


#### Models PREDICTS DATA ####
sub.sites.pred <- sites.pred # No further subsetting
sub.sites.pred$PREDICTS.LU <- droplevels(sub.sites.pred$PREDICTS.LU)
sub.sites.pred$logpop <- log1p(sub.sites.pred$pop)# logtransform population
sub.sites.pred$AvgTemp <- (sub.sites.pred$bio1_wc30s * .1)
sub.sites.pred$AvgPr <- (sub.sites.pred$bio12_wc30s)

# droplevels for Land use and land use-intensity
sub.sites.pred$PREDICTS.LUI[which(sub.sites.pred$PREDICTS.LUI=="Cannot decide")] <- NA
sub.sites.pred$PREDICTS.LU <- droplevels(sub.sites.pred$PREDICTS.LU)
sub.sites.pred$PREDICTS.LUI <- droplevels(sub.sites.pred$PREDICTS.LUI)
sub.sites.pred$UI <- factor(paste0(sub.sites.pred$PREDICTS.LU,".",sub.sites.pred$PREDICTS.LUI))
sub.sites.pred$LUInter <- (interaction(sub.sites.pred$PREDICTS.LU,sub.sites.pred$PREDICTS.LUI))
sub.sites.pred$LUInter <- droplevels(sub.sites.pred$LUInter)

sites$logpop <- log1p(sites$pop)
sites$PREDICTS.LU <- droplevels(sites$PREDICTS.LU)
sites$PREDICTS.LUI <- droplevels(sites$PREDICTS.LUI)
sites$AvgTemp <- (sites$bio1_wc30s * .1)
sites$AvgPr <- (sites$bio12_wc30s)
sites$Species_richness <- sites$spec
sites$LUInter <- interaction(sites$PREDICTS.LU,sites$PREDICTS.LUI)
sites$UI <- factor(paste0(sites$PREDICTS.LU,".",sites$PREDICTS.LUI))
sites$PREDICTS.LU <- droplevels(sites$PREDICTS.LU)
sites$PREDICTS.LUI <- droplevels(sites$PREDICTS.LUI)
sites$LUInter <- droplevels(sites$LUInter)


stop("Full stop if sourced!")

#### GO
response <- "logabund" # "Species_richness"
fam <- "gaussian" # "poisson"

p_full <- as.formula(paste(response," ~ ",
                          "PREDICTS.LU*PREDICTS.LUI +",
                          "Ecoregion +",
                          "elev +",
                          "meanEVI +",
                          "yieldEVI +",
                          "logpop + ",
                          "FC2000 +",
                          "AvgTemp +",
                          "AvgPr"
        )
)

# Landuse simple
p_lu <- as.formula(paste(response," ~ ",
                           "PREDICTS.LU*PREDICTS.LUI"                           
  )
)

# Only GIS predictors
p_gis <- as.formula(paste(response," ~ ",
                           "elev + ",
                           "logpop + ",
                           "meanEVI +",
                           "yieldEVI + ",
                           "FC2000"
                           ))

# Only Climate and topography
p_clim <- as.formula(paste(response," ~ ",
                           "AvgTemp +",
                           "AvgPr +",
                           "elev"
                           ))

  
# Predicts Nullmodel
p_null <- as.formula(paste(response," ~ ",
                           "1"
))

# Only modis derived  factors
p_modis <- as.formula(paste(response,"~",
                            "meanNDVI + ",
                            "meanEVI +",
                            "yieldNDVI +",
                            "yieldEVI"
                            )
                      )

#### LME 4 fitting ####
# calculate best random
rc1 <-  compare_randoms(sub.sites.pred,responseVar = response,fitFamily = fam,siteRandom = T,fixedFactors = c("Grouping"))

# Normal models
pfit <- glmer(update.formula(p_full,paste0("~ . +",rc1$best.random)),data=sub.sites.pred,family=fam,REML=F)#,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
pfit_null <- glmer(update.formula(p_null,paste0("~ . +",rc1$best.random)),data=sub.sites.pred,family=fam,REML=F) # fit nullmodel
pfit_lu <- glmer(update.formula(p_lu,paste0("~ . +",rc1$best.random)),data=sub.sites.pred,family=fam,REML=F) # Fit simple land-use model
pfit_gis <- glmer(update.formula(p_gis,paste0("~ . +",rc1$best.random)),data=sub.sites.pred,family=fam,REML=F)
pfit_clim <- glmer(update.formula(p_clim,paste0("~ . +",rc1$best.random)),data=sub.sites.pred,family=fam,REML=F)
pfit_mod <- glmer(update.formula(p_modis,paste0("~ . +",rc1$best.random)),data=sub.sites.pred,family=fam,REML=F)

aictab(cand.set=list(pfit,pfit_null,pfit_lu,pfit_gis,pfit_clim,pfit_mod),modnames=c("Full","Null","LU only","GIS metrics","Climate and Topo","MODIS"),second.ord=F,sort=T)
r.squaredLR(pfit,null=pfit_null) # Coeffient of determination compared to nullmodel
r.squaredGLMM(pfit)

# In case of non-convergance
gg <- pfit@optinfo$derivs$grad
hh <- pfit@optinfo$derivs$Hessian
vv <- sqrt(diag(solve(hh/2)))
print(summary(abs(gg*vv)))
#  the numbers printed here should be very small (examples that the developers deemed acceptable were around e-05.
relgrad <- with(pfit@optinfo$derivs,solve(Hessian,gradient))
print(max(abs(relgrad))) # this number should ideally be < 0.001.                                         
#source("AllFit.R")
#a <- allFit(pfit) # Try different optimizers
#summary.allfit(a)$which.OK

upa <- lmerTest::step(pfit,direction = "back")
upa <- glmer(upa$model@call,data = sub.sites.pred,family=fam,REML=F)
anova(pfit,upa,"chisq") # Sig differences
plot(Effect("meanEVI",upa))
pa <- aictab(list(upa,pfit_null),c("Updated model","Null-model"),second.ord = T,sort = T)
pa$r.squaredLR <- c(r.squaredLR(upa,pfit_null),0);pa <- as.data.frame(pa)
upa@call$formula

ups <- lmerTest::step(pfit,direction = "back")
ups <- glmer(ups$model@call,data = sub.sites.pred,family=fam,REML=F)
anova(pfit,ups,"chisq") # Sig differences
plot(allEffects(pfit_lu))
ps <- aictab(list(pfit_lu,pfit_null),c("Updated model","Null-model"),second.ord = T,sort = T)
ps$r.squaredLR <- c(r.squaredLR(pfit_lu,pfit_null),0);ps <- as.data.frame(ps)
ups@call$formula

# Export
pa <- pa[,c(3,4,9)];rownames(pa) <- NULL
ps <- ps[,c(3,4,9)];rownames(ps) <- NULL
names(pa) <- c("AICc","delta_AICc","r.squareLR")
names(ps) <- c("AICc","delta_AICc","r.squareLR")
rownames(pa) <- c("Land-use * Land-use intensity","Null-model")
rownames(ps) <- c("Land-use * Land-use intensity","Null-model")

library(gridExtra)
png("OutputComparison/ModelResultsPredicts1.png",width=8,height=6,units="in",type = "cairo",res=400)
p1 <- tableGrob(round(pa,3))
grid.arrange(p1)
dev.off()
png("OutputComparison/ModelResultsPredicts2.png",width=8,height=6,units="in",type = "cairo",res=400)
p2 <- tableGrob(round(ps,3))
grid.arrange(p2)
dev.off()

# --> Diagnostics <--
library(psych)
pairs.panels(sub.sites.pred[,c(37:46)],hist.col = "grey80")

mcp.fnc(pfit)
plot(pfit, resid(.)~fitted(.),type=c("p","smooth"),col="black",pch="+") # Fitted values against residuls, looks good
qqnorm(resid(pfit));qqline(resid(pfit),col="blue",lwd=2)
# Look for overdispersion
cbind(sum(residuals(pfit,type="pearson")^2),df.residual(pfit))
# Sum of squared Pearson residuals is way less than the residual degrees of freedom,
# thus not overdispersed

# Infestigate Effects
eff <- allEffects(pfit_clim)
plot(eff)

# Build successive model fitting
fit_GLMER <-function(data,response,predictorsNames,best_random,family,interactions=F,AICc=T){
  print("..Fitting sequential Mixed Effects models and select the best model...")
  # Save finished models
  ml <- list()
  ml_names <- vector()
  
  # Possible Combinations of Predictors
  for(id in 1:length(predictorsNames)){
    print(paste0("Adding ",id," predictors to models"))
    comb_mat <- combn(predictorsNames,id) # Build all combinations for x predictors
    for(i in 1:ncol(comb_mat)){
      b = comb_mat[,i]
      if( (id>1) & ("LUInter"%in%b)) b <- b[-which(b%in%c("PREDICTS.LU","PREDICTS.LUI"))]
      if(length(b)==0) break
      print(paste("Number of predictors in formula: ",length(b)))
      pred <- paste(b,collapse=" + ")  
      if(interactions&id>1) pred_int <- paste(b,collapse=" * ")      
      f <- as.formula(paste(response," ~ ",
                            pred," + ",
                            ifelse(interactions==T &id>1,paste0(pred_int," + "),""),
                            best_random # Add the best rando at the end
      ))
      print(f)
      try(fit <- glmer(f,data=data,family=family,REML=F))  
      ml <- c(ml,fit)
      ml_names <- c(ml_names,ifelse(interactions==T,paste0(pred," + ",pred_int),pred))
      print(paste0("Finished. Nr of Models: ",length(ml),".Next..."))
      rm(b)
    }
  }
  print("Make Model selection Table")
  o <- aictab(cand.set = ml,modnames = ml_names,second.ord = AICc,sort = T)
  o <- as.data.frame(o)
  #o$marg.rsquare <- unlist(lapply(ml,function(x) r.squaredGLMM(x)[1]))
  rm(ml,ml_names,id,fit,pred,pred_int)
  print("done!");print("-----------------------------------------")
  return(o)
}

fit_GLMER2 <-function(data,response,predictorsNames,best_random,family,AICc=T){
  print("..Fitting sequential Mixed Effects models and select the best model...")
  # Save finished models
  ml <- list()
  ml_names <- vector()
  
  # Possible Combinations of Predictors
  for(id in 1:length(predictorsNames)){
    print(paste0("Adding ",id," predictors to models"))
    comb_mat <- combn(predictorsNames,id) # Build all combinations for x predictors
    for(i in 1:ncol(comb_mat)){
      print(paste("Number of predictors in formula: ",length(comb_mat[,i])))
      pred <- paste(comb_mat[,i],collapse=" + ")  
      f <- as.formula(paste(response," ~ ",
                            "PREDICTS.LU*PREDICTS.LUI + ",
                            pred," + ",
                            best_random # Add the best rando at the end
      ))
      print(f)
      try(fit <- glmer(f,data=data,family=family,REML=F))  
      ml <- c(ml,fit)
      ml_names <- c(ml_names,pred)
      print(paste0("Finished. Nr of Models: ",length(ml),".Next..."))
    }
  }
  print("Make Model selection Table")
  o <- aictab.AICglmerMod(cand.set = ml,modnames = ml_names,second.ord = AICc,sort = T)
  o <- as.data.frame(o)
  #o$marg.rsquare <- unlist(lapply(ml,function(x) r.squaredGLMM(x)[1]))
  rm(ml,ml_names,id,fit,pred,pred_int)
  print("done!");print("-----------------------------------------")
  return(o)
}

# Kickout sites that are NA in PREDICTS.LU ?
#sub.sites.pred <- sub.sites.pred[which(!is.na(sub.sites.pred$PREDICTS.LU)),]


####HERE - model fits ####

# calculate best random
rc <-  CompareRandoms(sub.sites.pred,responseVar = "logabund",fitFamily = "gaussian",
                       fixedFactors = c("PREDICTS.LU"),siteRandom = F)
o <- fit_GLMER(sub.sites.pred,"logabund",c("PREDICTS.LU","PREDICTS.LUI","LUInter","logpop","yield.ndvi.corr","FC2000","meanNDVI"),rc$best.random,"gaussian",F,F)

rc2 <-  CompareRandoms(sub.sites.pred,responseVar = "Species_richness",fitFamily = "poisson",
                       fixedFactors = c("PREDICTS.LU"),siteRandom = T)
s <- fit_GLMER(sub.sites.pred,"Species_richness",c("PREDICTS.LU","PREDICTS.LUI","LUInter","logpop","yield.ndvi.corr","FC2000","meanNDVI"),rc2$best.random,"poisson",F,F)
# Maybe include Interaction manually

fit_a <- glmer(as.formula(paste("logabund ~ ",o$Modnames[1],"+",rc$best.random)),data = sub.sites.pred,family="gaussian",REML=F)
#sub.sites.pred$UI[which(grepl("NA",sub.sites.pred$UI))] <- NA;sub.sites.pred$UI <- droplevels(sub.sites.pred$UI)
#fit_a2 <- glmer("logabund~ UI + logpop + yield.ndvi.corr + meanNDVI + (1+PREDICTS.LU|SS)+(1|SSB)",data=sub.sites.pred,family="gaussian",REML=F)
fit_alu <- glmer(as.formula(paste("logabund ~ ","PREDICTS.LU","+",rc$best.random)),data = sub.sites.pred,family="gaussian",REML=F)
fit_s <- glmer(as.formula(paste("Species_richness ~",s$Modnames[1],"+",rc2$best.random)),data = sub.sites.pred,family="poisson",REML=F)
fit_slu <- glmer(as.formula(paste("Species_richness ~","PREDICTS.LU","+",rc2$best.random)),data = sub.sites.pred,family="poisson",REML=F)
fit_an <- glmer(as.formula(paste("logabund ~","1","+",rc$best.random)),data = sub.sites.pred,family="gaussian",REML=F)
fit_sn <- glmer(as.formula(paste("Species_richness ~","1","+",rc2$best.random)),data = sub.sites.pred,family="poisson",REML=F)
# Interaction
fit_alui <- glmer(as.formula(paste("logabund ~ ","LUInter","+",rc$best.random)),data = sub.sites.pred,family="gaussian",REML=F)
fit_slui <- glmer(as.formula(paste("Species_richness ~","LUInter","+",rc2$best.random)),data = sub.sites.pred,family="poisson",REML=F)

#AIC(fit_a,fit_a2)


vif.mer2 <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}
# Variance inflation factors
vif.mer2(fit_slu)

#m1<-glmer(Species_richness~PREDICTS.LU+PREDICTS.LUI+PREDICTS.LU:PREDICTS.LUI+(1|SS),data=sub.sites.pred,family="poisson")
#m2<-glmer(Species_richness~UI+(1|SS),data=sub.sites.pred,family="poisson") # drop NA
#m3<-glmer(Species_richness~interaction(PREDICTS.LU,PREDICTS.LUI)+(1|SS),data=sub.sites.pred,family="poisson")
#AIC(m1,m2,m3) ## all the same

r2.nakagawa.lme4(fit_a)$marginal - r2.nakagawa.lme4(fit_alu)$marginal

r2.nakagawa.lme4(fit_s)$marginal - r2.nakagawa.lme4(fit_slu)$marginal

# With interaction
r2.nakagawa.lme4(fit_a)$marginal - r2.nakagawa.lme4(fit_alui)$marginal
r2.nakagawa.lme4(fit_s)$marginal - r2.nakagawa.lme4(fit_slui)$marginal

a <- aictab(list(fit_a,fit_alu,fit_an),c(as.character(o$Modnames[1]),"Land use only","Null"),second.ord = F,sort=T);df.residual(fit_a);df.residual(fit_alu)
a <- aictab(list(fit_s,fit_slu,fit_sn),c(as.character(s$Modnames[1]),"Land use only","Null"),second.ord = F,sort=T);df.residual(fit_s);df.residual(fit_slu)

plot(allEffects(fit_s,se=T))
plot(Effect("logpop",fit_a))

# Evidence
anova(fit_alu,fit_a)
anova(fit_slu,fit_s)

# ---- #
## Test if slopes of models with Land-use only are different for broad-scale and independent ##
fit_alu <- glmer(as.formula(paste("logabund ~ ","PREDICTS.LU","+",rc$best.random)),data = sub.sites.pred,family="gaussian",REML=F)
fit_slu <- glmer(as.formula(paste("Species_richness ~","PREDICTS.LU","+",rc2$best.random)),data = sub.sites.pred,family="poisson",REML=F)
fit_a2 <- glmer(logabund ~ PREDICTS.LU + (1 | Transect),data = sub.sites,family="gaussian",REML=F)
fit_s2 <- glmer(Species_richness ~ PREDICTS.LU + (1 | Transect),data = sub.sites,family="poisson",REML=F)

## Abundance
# Z statistics for comparing slopes
# Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). Applied multiple regression/correlation analysis for the behavioral sciences (3rd ed.). Mahwah, New Jersey
b1 = summary(fit_alu)$coefficients[2]
b1_se = summary(fit_alu)$coefficients[7]
b2 = summary(fit_a2)$coefficients[2]
b2_se = summary(fit_a2)$coefficients[7]
( Z = (b1 - b2) / sqrt(( b1_se*b1_se + b2_se*b2_se )) )


## Species richness
b1 = summary(fit_slu)$coefficients[1]
b1_se = summary(fit_slu)$coefficients[6]
b2 = summary(fit_s2)$coefficients[1]
b2_se = summary(fit_s2)$coefficients[6]
Z = (b2 - b1) / sqrt(( b1_se*b1_se + b2_se*b2_se ))
abs(Z)

pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(modelFa)),mm))
fmpred_a <- data.frame(est=mm %*% fixef(modelFa),se=sqrt(pvar1))
fmpred_a$est.upp <- fmpred_a$est+sqrt(pvar1)#*1.96

# Z-statistics based on full plotting container
# Create first in customFigures
library(tidyr)
full2 <- full %>% filter(type == "log(Abundance)") %>% 
  mutate(se = abs(se.high - eff)) %>%  #recreate se
  dplyr::select(LandUse, eff,se, data) %>% # deselect
  group_by(LandUse,data) %>% 
  summarize(eff = mean(eff),
            se = mean(se)) # summarize
d1 = full2 %>% select(-se) %>% spread(data,value = eff) %>% rename(Field.eff = Field, All.eff = All) 
d2 = full2 %>% select(-eff) %>% spread(data,value = se) %>% rename(Field.se = Field, All.se = All) 
d <- left_join(d1,d2)  %>% 
  #Calculate Z score
  mutate(Z = (Field.eff - All.eff) / sqrt((All.se * All.se) + (Field.se * Field.se))  )

full2 <- full %>% filter(type == "Species Richness") %>% 
  mutate(se = abs(se.high - eff)) %>%  #recreate se
  dplyr::select(LandUse, eff,se, data) %>% # deselect
  group_by(LandUse,data) %>% 
  summarize(eff = mean(eff),
            se = mean(se)) # summarize
d1 = full2 %>% select(-se) %>% spread(data,value = eff) %>% rename(Field.eff = Field, All.eff = All) 
d2 = full2 %>% select(-eff) %>% spread(data,value = se) %>% rename(Field.se = Field, All.se = All) 
d <- left_join(d1,d2)  %>% 
  #Calculate Z score
  mutate(Z = (Field.eff - All.eff) / sqrt((All.se * All.se) + (Field.se * Field.se))  )


# ---- #
# Variable importance and evidence
# Sum up all aicwt of a balanced candidate model set that includes this model
# Estimates of the relative importance of predictor variables xj can best be made by
# summing the AIC weights across all the models in the set where variable j occurs.
# The larger the sum of model weights, the more important a variable is relative to the other variables. 
# The relative importance values can be used to rank the variables. 
# However, to use this method, one must have an equal number of models for each variable; 
# otherwise, some variables will be over represented or under represented resulting in biased relative importance values.

# http://stats.stackexchange.com/questions/25322/relative-variable-importance-values-vs-magnitude-of-effect

# abundance
sum(o$AICWt[grep("logpop",o$Modnames)])
sum(o$AICWt[grep("FC2000",o$Modnames)])
sum(o$AICWt[grep("meanNDVI",o$Modnames)])
sum(o$AICWt[grep("yield.ndvi.corr",o$Modnames)])
# LU
sum(o$AICWt[grep("LUInter",o$Modnames,invert = T)&& grep("PREDICTS.LU",o$Modnames)]) # all LU without LUInter
sum(o$AICWt[grep("LUInter",o$Modnames,invert = T)&& grep("PREDICTS.LUI",o$Modnames)])
sum(o$AICWt[grep("LUInter",o$Modnames)&& grep("PREDICTS.LU",o$Modnames,invert = T)&& grep("PREDICTS.LUI",o$Modnames,invert = T)])

# Species richness
sum(s$AICWt[grep("logpop",s$Modnames)])
sum(s$AICWt[grep("FC2000",s$Modnames)])
sum(s$AICWt[grep("meanNDVI",s$Modnames)])
sum(s$AICWt[grep("yield.ndvi.corr",s$Modnames)])
# Mod
sum(s$AICWt[grep("LUInter",s$Modnames,invert = T)&& grep("PREDICTS.LUI",s$Modnames)]) # all LU without LUInter
sum(s$AICWt[grep("LUInter",s$Modnames,invert = T)&& grep("PREDICTS.LUI",s$Modnames)]) # all LUI without LUInter
sum(s$AICWt[grep("LUInter",s$Modnames)&& grep("PREDICTS.LU",s$Modnames,invert = T)&& grep("PREDICTS.LUI",s$Modnames,invert = T)])


#### Spatial autocorrelation - Histos####
# Make a histogram for abundance and richness
# like in Newbold 2015 et al. and save as appendix
SpaTest <- function(model,all.data){
  model.data <- model@frame
  model.data$res <- residuals(model)
  if ("SSBS"%in% names(model.data)){
    model.data$Longitude <- all.data$Longitude[match(model.data$SSBS, 
                                                     all.data$SSBS)]
    model.data$Latitude <- all.data$Latitude[match(model.data$SSBS, 
                                                   all.data$SSBS)]
  } else if("SSB"%in% names(model.data)){
    model.data$Longitude <- all.data$Longitude[match(model.data$SSB, 
                                                     all.data$SSB)]
    model.data$Latitude <- all.data$Latitude[match(model.data$SSB, 
                                                   all.data$SSB)]
  }
  studies <- character()
  failed <- character()
  moran.i <- numeric()
  moran.p <- numeric()
  i = 1
  for (ss in unique(model.data$SS)) {
    cat(paste("\rProcessing study ", i, " of ", length(unique(model.data$SS)), 
              sep = ""))
    data.sub <- model.data[model.data$SS == ss, ]
    ds.nb <- try(dnearneigh(cbind(data.sub$Longitude, data.sub$Latitude), 
                            d1 = 1e-08, d2 = 10), silent = TRUE)
    ds.listw <- try(nb2listw(ds.nb), silent = TRUE)
    mt <- try(moran.test(data.sub$res, ds.listw), silent = TRUE)
    if (class(mt) == "htest") {
      if ((!is.na(mt$statistic))) {
        studies <- c(studies, ss)
        moran.i <- c(moran.i, mt$statistic)
        moran.p <- c(moran.p, mt$p.value)
      }
      else {
        failed <- c(failed, ss)
      }
    }
    else {
      failed <- c(failed, ss)
    }
    i <- i + 1
  }
  return(list(studies = studies, I = moran.i, P = moran.p, 
              failed = failed))
}

spar_a <- SpaTest(fit_a,sub.sites.pred) # MJ1_2013__Ndanganga 2
spar_s <- SpaTest(fit_s,sub.sites.pred)

d <- data.frame(Type=character(),P=numeric(),I=numeric())
d <- rbind(d,data.frame(Type="Abundance",P=spar_a$P,I=spar_a$I))
d <- rbind(d,data.frame(Type="Species Richness",P=spar_s$P,I=spar_s$I))
d <- subset(d,!is.infinite(d$I)) # Remove infs
d$AC <- d$P <=0.05
# Plot
library(ggthemr)
ggthemr("greyscale",layout="clean",type="outer",spacing=0)
g <- ggplot(data=d,aes(x=P),fill="grey")
g <- g + geom_histogram() + geom_vline(aes(xintercept=0.05),color="red")
g <- g + facet_wrap(~Type,as.table = T,nrow = 2)
g <- g + scale_x_continuous(limits=c(0,1),breaks=pretty_breaks(n=10))
g <- g + labs(x = "P-value", y= "Frequency")
g
ggsave("../000Thesis_writeup/PaperDraft/Appendix/SI Figure2.png",plot=g,scale = 1,units = "mm",dpi=400)

spar <- SpaTest(fit_alu,sub.sites.pred) # MJ1_2013__Ndanganga 2
spar <- SpaTest(fit_slu,sub.sites.pred) # "MJ1_2008__Munyekenye 1"

# RMSE of abundance and SR
re <- data.frame(PREDICTS.LU=sub.sites.pred$PREDICTS.LU,PREDICTS.LUI=sub.sites.pred$PREDICTS.LUI)
re$resid_a <- resid(fit_a)[match(row.names(sub.sites.pred),names(resid(fit_a)))]
re$resid_s <- resid(fit_s)[match(row.names(sub.sites.pred),names(resid(fit_s)))]

ddply(re,.(PREDICTS.LU,PREDICTS.LUI),summarise,
      RMSE_a=round( sqrt(mean(resid_a,na.rm = T)^2),3 ),
      RMSE_s=round( sqrt(mean(resid_s,na.rm = T)^2),3 )
      )

#### LME 4 Prediction and Results ####
# display results of models
#apply to field sties
#log-transform predicted values, calculate average and backtransform

#coefficents apart from reference level

#if predictors don't improve fit, other env. are needed
#or traits are different.

# Figure 5 Relative Divide intercept through ###
sitemakeRatios2 <- function(sub.sites,measure,pred,exp=T,Grouping="Birds") {
  fw_ab <- summarySEwithin(sub.sites,measurevar=measure,withinvars=c(pred))
  fw_ab$ind <- match(fw_ab[,1],ord)
  fw_ab$se.high <- fw_ab$val + fw_ab$se  
  fw_ab$se.low <- fw_ab$val - fw_ab$se
  if(exp){
    fw_ab$val <- exp(fw_ab$val)
    fw_ab$se.high <- exp(fw_ab$se.high)
    fw_ab$se.low <- exp(fw_ab$se.low)
  }
  # Reference
  t <- fw_ab$val[which(fw_ab$ind==1)]
  
  fw_ab$nratio <- (fw_ab$val/t )*100
  
  fw_ab$se.low.n <- (fw_ab$se.low /fw_ab$se.low[which(fw_ab$ind==1)] )*100
  fw_ab$se.high.n <- (fw_ab$se.high /fw_ab$se.high[which(fw_ab$ind==1)])*100
  
  of <- data.frame(LandUse=fw_ab[,1],eff=fw_ab$nratio,se.low=fw_ab$se.low.n,se.high=fw_ab$se.high.n,data="Field",Grouping=Grouping)
  of$LandUse <- droplevels(ordered(of$LandUse,levels=(ord)))
  return(of)
}

# Take Sites data predictions
# Kick urban intense out for now because it is not present in the broad-scale dataset
sa <- sites[which(sites$LUInter!="Urban.Intense use"),]
sa$LUInter <- droplevels(sa$LUInter)
sa$PREDICTS.LU <- droplevels(sa$PREDICTS.LU);sa$PREDICTS.LUI <- droplevels(sa$PREDICTS.LUI)
sa <- subset(sa,Grouping=="Aves")  # Kick out Dickens
# Insert some dummy variables for the random intercepts
sa$Species_richness <- sa$spec
sa$SS <- sa$Transect;sa$SSB <- paste0(sa$SS,"_1_");sa$SSBS <- paste0(sa$SSB,sa$ID)

#take best model, apply predict absolute abundance
w <- names(fit_a@frame)#;w <- w[-which(w%in%c("logabund"))]
grd_abd <- subset(sa,select=w)

#RmNonObserved <- function(fit,pred,fam="gaussian"){
#to <- unique(fit@frame$LUInter)[which(!unique(fit@frame$LUInter) %in% (unique(pred$LUInter)))]
#if(length(to)>0){
#  # Set LUInter to NA for not observed interactions and refit
#  altp <- sub.sites.pred
#  altp$LUInter[which(altp$LUInter%in%to)] <- NA
#  altp$LUInter <- droplevels(altp$LUInter)
#  fit <- glmer(fit@call$formula,data = altp,family=fam,REML=F)
#  RmNonObserved(fit,pred)
#} else return(fit)
#}
#fit_a <- RmNonObserved(fit_a,grd_abd)

# Find the levels
all.levels.PREDICTS.LU<-union(levels(grd_abd$PREDICTS.LU),levels(fit_a@frame$PREDICTS.LU))
all.levels.PREDICTS.LUI<-union(levels(grd_abd$PREDICTS.LUI),levels(fit_a@frame$PREDICTS.LUI))
all.levels.LUInter<-union(levels(grd_abd$LUInter),levels(fit_a@frame$LUInter))
# Make modelling data where the factors have the same levels
model.data<-fit_a@frame
model.data$PREDICTS.LU<-factor(model.data$PREDICTS.LU,levels=all.levels.PREDICTS.LU)
model.data$LUInter <- factor(model.data$LUInter,levels=all.levels.LUInter)
#model.data$PREDICTS.LUI<-factor(model.data$PREDICTS.LUI,levels=all.levels.PREDICTS.LUI)
# Rebuild the model with the new data
modelFa<-glmer(formula = as.formula(paste("logabund ~",o$Modnames[1],"+",rc$best.random)), data = model.data, family = "gaussian")
modelLa<-glmer(as.formula(paste("logabund ~ ","PREDICTS.LU","+",rc$best.random)),data = model.data,family="gaussian",REML=F)
# Make prediction data with the corrected levels
grd_abd$PREDICTS.LU<-factor(grd_abd$PREDICTS.LU,levels=all.levels.PREDICTS.LU)
#grd_abd$PREDICTS.LUI<-factor(grd_abd$PREDICTS.LUI,levels=all.levels.PREDICTS.LUI)
grd_abd$LUInter <- factor(grd_abd$LUInter,levels=all.levels.LUInter)
# Now predict with the model
mm<-model.matrix(terms(modelFa),grd_abd)
pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(modelFa)),mm))
fmpred_a <- data.frame(est=mm %*% fixef(modelFa),se=sqrt(pvar1))
fmpred_a$est.upp <- fmpred_a$est+sqrt(pvar1)#*1.96
fmpred_a$est.low <- fmpred_a$est-sqrt(pvar1)#*1.96
mm<-model.matrix(terms(modelLa),grd_abd)
pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(modelLa)),mm))
fmpred_alu <- data.frame(est=mm %*% fixef(modelLa),se=sqrt(pvar1))
fmpred_alu$est.upp <- fmpred_alu$est+sqrt(pvar1)#*1.96
fmpred_alu$est.low <- fmpred_alu$est-sqrt(pvar1)#*1.96

#fpred_a <- predict(modelFa,newdata=grd_abd,type="link",re.form=NA,allow.new.levels=F)#,re.form=NA)
#fpred_alu <- predict(modelLa,newdata=grd_abd,type="link",re.form=NA,allow.new.levels=F)#,re.form=NA)

pd1 <- data.frame(id=sa$ID,LU=sa$PREDICTS.LU,val=fmpred_a$est,set=fmpred_a$se,val.se.upp=fmpred_a$est.upp,val.se.low=fmpred_a$est.low,Type="(log) Abundance",Transect=sa$Transect,Grouping=sa$Grouping)
pd2 <- data.frame(id=sa$ID,LU=sa$PREDICTS.LU,val=fmpred_alu$est,set=fmpred_alu$se,val.se.upp=fmpred_alu$est.upp,val.se.low=fmpred_alu$est.low,Type="(log) Abundance",Transect=sa$Transect,Grouping=sa$Grouping)

# Full Model
detachPackage("dplyr")
fw_ab <- summarySEwithin(pd1,measurevar="val",withinvars=c("LU","Transect"))
fw_ab$ind <- match(fw_ab[,1],ord)
fw_ab$se <- summarySEwithin(pd1,measurevar="set",withinvars=c("LU","Transect"))$set # replace with aggregated prediction error
fw_ab$se.high <- fw_ab$val + fw_ab$se
fw_ab$se.low <- fw_ab$val - fw_ab$se
fw_ab$se.high.n <-NA;fw_ab$se.low.n <- NA
# Reference
t <- fw_ab$val[which(fw_ab$ind==1&fw_ab$Transect=="Taita")]
fw_ab$nratio[which(fw_ab$Transect=="Taita")] <- (fw_ab$val[which(fw_ab$Transect=="Taita")] /t )*100
fw_ab$se.low.n[which(fw_ab$Transect=="Taita")] <- (fw_ab$se.low[which(fw_ab$Transect=="Taita")] /fw_ab$se.low[which(fw_ab$ind==1&fw_ab$Transect=="Taita")] )*100
fw_ab$se.high.n[which(fw_ab$Transect=="Taita")] <- (fw_ab$se.high[which(fw_ab$Transect=="Taita")] /fw_ab$se.high[which(fw_ab$ind==1&fw_ab$Transect=="Taita")])*100

k <- fw_ab$val[which(fw_ab$ind==1&fw_ab$Transect=="Kilimanjaro")]
fw_ab$nratio[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$val[which(fw_ab$Transect=="Kilimanjaro")] /k )*100
fw_ab$se.low.n[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$se.low[which(fw_ab$Transect=="Kilimanjaro")] /fw_ab$se.low[which(fw_ab$ind==1&fw_ab$Transect=="Kilimanjaro")] )*100
fw_ab$se.high.n[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$se.high[which(fw_ab$Transect=="Kilimanjaro")] /fw_ab$se.high[which(fw_ab$ind==1&fw_ab$Transect=="Kilimanjaro")])*100
dfc <- data.frame(LandUse=fw_ab[,1],Transect=fw_ab$Transect,eff=fw_ab$nratio,se.low=fw_ab$se.low.n,se.high=fw_ab$se.high.n)
dfc$LandUse <- droplevels(ordered(dfc$LandUse,levels=(ord)))
dfc$data<-"Best Model";dfc$type <- "(log) Abundance"
rm(t,k,fw_ab)

# LandUse only
fw_ab <- summarySEwithin(pd2,measurevar="val",withinvars=c("LU","Transect"))
fw_ab$ind <- match(fw_ab[,1],ord)
fw_ab$se <- summarySEwithin(pd2,measurevar="set",withinvars=c("LU","Transect"))$set # replace with aggregated prediction error
fw_ab$se.high <- fw_ab$val + fw_ab$se
fw_ab$se.low <- fw_ab$val - fw_ab$se
fw_ab$se.high.n <-NA;fw_ab$se.low.n <- NA
# Reference
t <- fw_ab$val[which(fw_ab$ind==1&fw_ab$Transect=="Taita")]
fw_ab$nratio[which(fw_ab$Transect=="Taita")] <- (fw_ab$val[which(fw_ab$Transect=="Taita")] /t )*100
fw_ab$se.low.n[which(fw_ab$Transect=="Taita")] <- (fw_ab$se.low[which(fw_ab$Transect=="Taita")] /fw_ab$se.low[which(fw_ab$ind==1&fw_ab$Transect=="Taita")] )*100
fw_ab$se.high.n[which(fw_ab$Transect=="Taita")] <- (fw_ab$se.high[which(fw_ab$Transect=="Taita")] /fw_ab$se.high[which(fw_ab$ind==1&fw_ab$Transect=="Taita")])*100

k <- fw_ab$val[which(fw_ab$ind==1&fw_ab$Transect=="Kilimanjaro")]
fw_ab$nratio[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$val[which(fw_ab$Transect=="Kilimanjaro")] /k )*100
fw_ab$se.low.n[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$se.low[which(fw_ab$Transect=="Kilimanjaro")] /fw_ab$se.low[which(fw_ab$ind==1&fw_ab$Transect=="Kilimanjaro")] )*100
fw_ab$se.high.n[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$se.high[which(fw_ab$Transect=="Kilimanjaro")] /fw_ab$se.high[which(fw_ab$ind==1&fw_ab$Transect=="Kilimanjaro")])*100
dfc2 <- data.frame(LandUse=fw_ab[,1],Transect=fw_ab$Transect,eff=fw_ab$nratio,se.low=fw_ab$se.low.n,se.high=fw_ab$se.high.n)
dfc2$LandUse <- droplevels(ordered(dfc2$LandUse,levels=(ord)))
dfc2$data<-"Land use only";dfc2$type <- "(log) Abundance"
rm(t,k,fw_ab)

d1 <- rbind(dfc,dfc2)

# species richness
w <- names(fit_s@frame)#;w <- w[-which(w%in%c("Species_richness","SS","SSB","SSBS"))]
grd_spe <- subset(sa,select=w)
#fit_s <- RmNonObserved(fit_s,grd_spe,"poisson")

# Find the levels
all.levels.PREDICTS.LU<-union(levels(grd_spe$PREDICTS.LU),levels(fit_s@frame$PREDICTS.LU))
#all.levels.PREDICTS.LUI<-union(levels(grd_spe$PREDICTS.LUI),levels(fit_s@frame$PREDICTS.LUI))
all.levels.LUInter<-union(levels(grd_spe$LUInter),levels(fit_s@frame$LUInter))
# Make modelling data where the factors have the same levels
model.data<-fit_s@frame
model.data$PREDICTS.LU<-factor(model.data$PREDICTS.LU,levels=all.levels.PREDICTS.LU)
#model.data$PREDICTS.LUI<-factor(model.data$PREDICTS.LUI,levels=all.levels.PREDICTS.LUI)
model.data$LUInter <- factor(model.data$LUInter,levels=all.levels.LUInter)
# Rebuild the model with the new data
modelFs<-glmer(formula = as.formula(paste("Species_richness ~",s$Modnames[1],"+",rc2$best.random)), data = model.data, family = "poisson")
modelLs<-glmer(as.formula(paste("Species_richness ~ ","PREDICTS.LU","+",rc2$best.random)),data = model.data,family="poisson",REML=F)
# Make prediction data with the corrected levels
grd_spe$PREDICTS.LU<-factor(grd_spe$PREDICTS.LU,levels=all.levels.PREDICTS.LU)
#grd_spe$PREDICTS.LUI<-factor(grd_spe$PREDICTS.LUI,levels=all.levels.PREDICTS.LUI)
grd_spe$LUInter <- factor(grd_spe$LUInter,levels=all.levels.LUInter)

# Now predict with the model
mm<-model.matrix(terms(modelFs),grd_spe)
pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(modelFs)),mm))
fmpred_s <- data.frame(est=mm %*% fixef(modelFs),se=sqrt(pvar1))
fmpred_s$est.upp <- fmpred_s$est+sqrt(pvar1)#*1.96
fmpred_s$est.low <- fmpred_s$est-sqrt(pvar1)#*1.96
mm<-model.matrix(terms(modelLs),grd_spe)
pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(modelLs)),mm))
fmpred_slu <- data.frame(est=mm %*% fixef(modelLs),se=sqrt(pvar1))
fmpred_slu$est.upp <- fmpred_slu$est+sqrt(pvar1)#*1.96
fmpred_slu$est.low <- fmpred_slu$est-sqrt(pvar1)#*1.96

#fpred_s <- predict(modelFs,newdata=grd_spe,type="link",re.form=NA,allow.new.levels=F)
#fpred_slu <- predict(modelLs,newdata=grd_spe,type="link",re.form=NA,allow.new.levels=F)

pds <- data.frame(id=sa$ID,LU=sa$PREDICTS.LU,set=fmpred_s$se,val=fmpred_s$est,val.se.upp=fmpred_s$est.upp,val.se.low=fmpred_s$est.low,Type="Species richness",Transect=sa$Transect,Grouping=sa$Grouping)
pds2 <- data.frame(id=sa$ID,LU=sa$PREDICTS.LU,set=fmpred_slu$se,val=fmpred_slu$est,val.se.upp=fmpred_slu$est.upp,val.se.low=fmpred_slu$est.low,Type="Species richness",Transect=sa$Transect,Grouping=sa$Grouping)

# Full Model
fw_ab <- summarySEwithin(pds,measurevar="val",withinvars=c("LU","Transect"))
fw_ab$ind <- match(fw_ab[,1],ord)
fw_ab$val <- exp(fw_ab$val)
fw_ab$se <- exp(summarySEwithin(pds,measurevar="set",withinvars=c("LU","Transect"))$set) # replace with aggregated prediction error
fw_ab$se.high <- fw_ab$val + fw_ab$se
fw_ab$se.low <- fw_ab$val - fw_ab$se
fw_ab$se.high.n <-NA;fw_ab$se.low.n <- NA;fw_ab$nratio <- NA
# Reference
t <- fw_ab$val[which(fw_ab$ind==1&fw_ab$Transect=="Taita")]
fw_ab$nratio[which(fw_ab$Transect=="Taita")] <- (fw_ab$val[which(fw_ab$Transect=="Taita")] /t )*100
fw_ab$se.low.n[which(fw_ab$Transect=="Taita")] <- (fw_ab$se.low[which(fw_ab$Transect=="Taita")] /fw_ab$se.low[which(fw_ab$ind==1&fw_ab$Transect=="Taita")] )*100
fw_ab$se.high.n[which(fw_ab$Transect=="Taita")] <- (fw_ab$se.high[which(fw_ab$Transect=="Taita")] /fw_ab$se.high[which(fw_ab$ind==1&fw_ab$Transect=="Taita")])*100

k <- fw_ab$val[which(fw_ab$ind==1&fw_ab$Transect=="Kilimanjaro")]
fw_ab$nratio[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$val[which(fw_ab$Transect=="Kilimanjaro")] /k )*100
fw_ab$se.low.n[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$se.low[which(fw_ab$Transect=="Kilimanjaro")] /fw_ab$se.low[which(fw_ab$ind==1&fw_ab$Transect=="Kilimanjaro")] )*100
fw_ab$se.high.n[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$se.high[which(fw_ab$Transect=="Kilimanjaro")] /fw_ab$se.high[which(fw_ab$ind==1&fw_ab$Transect=="Kilimanjaro")])*100
dfr <- data.frame(LandUse=fw_ab[,1],Transect=fw_ab$Transect,eff=fw_ab$nratio,se.low=fw_ab$se.low.n,se.high=fw_ab$se.high.n)
dfr$LandUse <- droplevels(ordered(dfr$LandUse,levels=(ord)))
dfr$data<-"Best Model";dfr$type<-"Species richness"
rm(t,k,fw_ab)

# LandUse only
fw_ab <- summarySEwithin(pds2,measurevar="val",withinvars=c("LU","Transect"))
fw_ab$ind <- match(fw_ab[,1],ord)
fw_ab$val <- exp(fw_ab$val)
fw_ab$se <- exp(summarySEwithin(pds2,measurevar="set",withinvars=c("LU","Transect"))$set) # replace with aggregated prediction error
fw_ab$se.high <- fw_ab$val + fw_ab$se
fw_ab$se.low <- fw_ab$val - fw_ab$se
fw_ab$se.high.n <-NA;fw_ab$se.low.n <- NA;fw_ab$nratio <- NA
# Reference
t <- fw_ab$val[which(fw_ab$ind==1&fw_ab$Transect=="Taita")]
fw_ab$nratio[which(fw_ab$Transect=="Taita")] <- (fw_ab$val[which(fw_ab$Transect=="Taita")] /t )*100
fw_ab$se.low.n[which(fw_ab$Transect=="Taita")] <- (fw_ab$se.low[which(fw_ab$Transect=="Taita")] /fw_ab$se.low[which(fw_ab$ind==1&fw_ab$Transect=="Taita")] )*100
fw_ab$se.high.n[which(fw_ab$Transect=="Taita")] <- (fw_ab$se.high[which(fw_ab$Transect=="Taita")] /fw_ab$se.high[which(fw_ab$ind==1&fw_ab$Transect=="Taita")])*100

k <- fw_ab$val[which(fw_ab$ind==1&fw_ab$Transect=="Kilimanjaro")]
fw_ab$nratio[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$val[which(fw_ab$Transect=="Kilimanjaro")] /k )*100
fw_ab$se.low.n[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$se.low[which(fw_ab$Transect=="Kilimanjaro")] /fw_ab$se.low[which(fw_ab$ind==1&fw_ab$Transect=="Kilimanjaro")] )*100
fw_ab$se.high.n[which(fw_ab$Transect=="Kilimanjaro")] <- (fw_ab$se.high[which(fw_ab$Transect=="Kilimanjaro")] /fw_ab$se.high[which(fw_ab$ind==1&fw_ab$Transect=="Kilimanjaro")])*100
dfr2 <- data.frame(LandUse=fw_ab[,1],Transect=fw_ab$Transect,eff=fw_ab$nratio,se.low=fw_ab$se.low.n,se.high=fw_ab$se.high.n)
dfr2$LandUse <- droplevels(ordered(dfr2$LandUse,levels=(ord)))
dfr2$data<-"Land use only";dfr2$type<-"Species richness"
rm(t,k,fw_ab)

d2 <- rbind(dfr,dfr2)

d <- rbind(d1,d2)
d$LandUse <- ordered(d$LandUse,ord)

rect_lab <- c("PV","SV","PL","CL","UR")
xmin = c(0.5,1.5,2.5,3.5,4.5);xmax=c(1.5,2.5,3.5,4.5,5.5)
# Make a plot
library(ggthemes)
g <- ggplot(d,aes(x=LandUse,y=eff,fill=data,shape=data,group=data)) + theme_light(base_size = 16)
g <- g + geom_hline(yintercept=100, lty=2,size=.5)
g <- g + geom_vline(xintercept=xmax[-5],size=.2,linetype="dotted") + geom_vline(xintercept=0)# and vline seperation
g <- g + geom_pointrange(aes(ymax=se.high,ymin=se.low),size=.9,position=position_jitterdodge(jitter.width = 1,jitter.height = 1,dodge.width = 0.5) )
g <- g + facet_wrap(Transect~type,scales = "free_y",as.table = F,drop = T) #+ theme(strip.text.x = element_text(size=16), strip.background = element_rect(colour="black"))
g <- g + scale_y_continuous(expand=c(0,1),breaks=pretty_breaks()) + coord_cartesian()
#g <- g + theme(axis.text.y=element_text(face = c(rep("plain",4),"bold",rep("plain",10)) )) # Special Eyecandy of 100% 
g <- g + scale_x_discrete(labels=rect_lab)
g <- g + labs(y="Percentage of Baseline",x="")
g <- g + ggtitle("Broad-scale model prediction \n on field study sites")
g <- g + scale_fill_discrete(guide="none")
g <- g + scale_shape_manual(values=c(15,16),labels = c("Best Model","Land use only"),
                            guide = guide_legend(title = "Dataset:",direction = "vertical", title.position = "left",
                                                 label.position="right", label.hjust = 0.5, label.vjust = 0.5,nrow=1)) 
g <- g + theme(legend.position="bottom",axis.title.x=element_blank(),axis.ticks.x=element_blank())
g
ggsave("OutputComparison/PredictionComparison.png",plot=g,dpi=400)

# Calculate the Prediction deviation from the observed for both models
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

sa$val <- sa$logabund
ofA <- sitemakeRatios(sa,"PREDICTS.LU",exp=T,"Birds")
ofA$type <- "(log) Abundance";ofA$Grouping <- "Field"

sa$val <- sa$spec
ofA_r <- sitemakeRatios(sa,"PREDICTS.LU",exp=F,"Birds")
ofA_r$type <- "Species richness";ofA_r$Grouping <- "Field"

d$Grouping <- "PREDICTS"
full <- rbind(ofA,ofA_r,d)
full$LandUse <- ordered(full$LandUse,ord)
full$Transect <- ordered(full$Transect,levels=c("Kilimanjaro","Taita"))
full$se <- abs(full$eff - full$se.low)
d2 <- ddply(full,.(Transect,type,LandUse),summarise,
      # Also format to Percentage scale 0-1
      Top= 0.01 * (eff[data=="Field"] - eff[data=="Best Model"]),
      Top.se = 0.01 * (se[data=="Field"] - se[data=="Best Model"]),
      LU= 0.01 * (eff[data=="Field"] - eff[data=="Land use only"]),
      LU.se = 0.01 * (se[data=="Field"] - se[data=="Land use only"])
      )
d2_se <- d2 %>% dplyr::select(Transect:LandUse,Top.se,LU.se) %>% melt() %>% mutate(variable = str_replace(variable,"\\.se",""))
d2_se<- d2_se[which(d2_se$LandUse!="Primary Vegetation"),] # Kickout Primary Vegetation

d2 <- d2 %>% dplyr::select(Transect:LandUse,Top,LU) %>%melt()
d2<- d2[which(d2$LandUse!="Primary Vegetation"),] # Kickout Primary Vegetation
# Join back
d2 <- d2_se %>% dplyr::rename(se = value) %>% left_join(d2,.)
d2$LandUse <- ordered(d2$LandUse,rev(ord))

# Z statistics for comparing slopes
# Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). Applied multiple regression/correlation analysis for the behavioral sciences (3rd ed.). Mahwah, New Jersey
z <-full %>% dplyr::filter(Grouping == "PREDICTS") %>% dplyr::select(LandUse,Transect, type,data, eff,se.low) %>% 
  mutate(se = abs(eff - se.low) ) %>% dplyr::select(-se.low) %>%
  dplyr::group_by(Transect,type,LandUse) %>% 
  summarise( 
        Z = (eff[data=="Best Model"] - eff[data=="Land use only"]) / sqrt(( se[data=="Best Model"]*se[data=="Best Model"] + se[data=="Land use only"]*se[data=="Land use only"] )) 
        ) 
fz <- full %>% dplyr::filter(Grouping == "PREDICTS") %>% dplyr::select(LandUse,Transect, type,data, eff,se.low) %>% 
  mutate(se = abs(eff - se.low) ) %>% dplyr::select(-se.low) %>%
  dplyr::select(LandUse:type) %>% dplyr::distinct()
fz$Z <- z

# Alternative:
# Calculate average difference per site,model and type
d2 %>% group_by(type,variable) %>% 
  summarise(avg = (mean(abs(value)))*100)

p <- ggplot(d2,aes(x=LandUse,y=value,group=variable,fill=factor(variable))) + theme_few(base_size = 16)
p <- p + geom_hline(yintercept=0)
p <- p + geom_bar(stat = "identity",position="dodge")
p <- p + geom_errorbar(aes(x = LandUse, ymin = value-se,ymax = value+se,group=variable),inherit.aes = F, position=position_dodge(.9),stat = "identity",width = 0.2)
p <- p + coord_flip()
p <- p + facet_grid(type~Transect,drop = T,scales = "free")
p <- p + labs(x="",y="Average difference from observed\n for both study transects")
p <- p + scale_y_continuous(breaks=pretty_breaks(5),labels=percent_format())
p <- p + theme(strip.text.x = element_text(size=11))#, strip.background = element_rect(colour="black"))
p <- p + scale_fill_manual(values = c("#3F3F3F","#BFBFBF"),labels = c("Land use only","Best selected"),
                            guide = guide_legend(title = "Model",direction = "vertical", title.position = "top",
                                                 label.position="right")) 
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1,face = "bold"))
p
ggsave("Figure4.png",plot=p,scale=1.1,dpi=400)
# For THESIS
#ggsave("../000Thesis_writeup/WriteUpLatex/gfx/Figure4.pdf",plot=p,scale=1.1,dpi=400)

# Transparent
t <- theme(plot.background = element_rect(fill = "transparent",colour = NA),panel.background = element_rect(fill = "transparent",colour = NA),
           axis.text = element_text(colour="white"),axis.text.x =element_text(colour="white"),axis.text.y =element_text(colour="white"),axis.ticks = element_line(colour = "white"), axis.line = element_line(colour="white"),axis.title = element_text(colour="white"),
           legend.key = element_rect(fill=NA,colour = NA),legend.background = element_rect(fill = NA,colour=NA),legend.text = element_text(colour="white"),
           title = element_text(colour="white"),
           strip.background = element_rect(fill=NA), strip.text = element_text(color="white"),
           panel.grid.major = element_blank(),panel.grid.minor = element_blank()
) 


ggsave(filename="../001Thesis_Defense/assets/img/Figure4_trans.png",plot=p+t,scale=1.1,dpi=400,bg = "transparent")


#### Hypothesis - Overperforming Model Comparison ####
# What makes overperforming studies all over Africa perform better.
# underpforming worse
# Are those factors the same as in my fieldwork case
# Idea -> get over and under outliers for the PREDICTS data (cooks distances)
sub.sites <- sites#subset(sites,PREDICTS.LU!="Urban") # Kick out urban for now
sub.sites$PREDICTS.LU <- droplevels(sub.sites$PREDICTS.LU)
sub.sites$logpop <- log(sub.sites$pop) # Logtransform Population
# Start constructing glm on field data using available predictors
library(psych)
pairs.panels(sub.sites[,c(55:60)],hist.col = "grey80")
ga_modis <- glm(spec~meanNDVI+meanEVI+yieldNDVI+yieldEVI,data=sub.sites,family=poisson)
up <- stepAIC(ga_modis,direction="both",trace=T) # all kicked out...

names(sub.sites)
ga_full <- glm(logabund~PREDICTS.LU*PREDICTS.LUI+Elev+I(Elev^2)+logpop+I(logpop^2)+yieldEVI+I(yieldEVI^2)+Transect,data=sub.sites,family=gaussian)
ga_lu   <- glm(logabund~PREDICTS.LU*PREDICTS.LUI+Transect,data=sub.sites,family=gaussian)
ga_rs   <- glm(logabund~Elev+I(Elev^2)+logpop+I(logpop^2)+yieldEVI+I(yieldEVI^2)+Transect,data=sub.sites,family=gaussian) 
ga_null <- glm(logabund~1,data=sub.sites,family = gaussian)
gs_full <- glm(spec~PREDICTS.LU*PREDICTS.LUI+Elev+I(Elev^2)+logpop+I(logpop^2)+yieldEVI+I(yieldEVI^2)+Transect,data=sub.sites,family=poisson(link=log)) # nothing to be dropped
gs_lu   <- glm(spec~PREDICTS.LU*PREDICTS.LUI+Transect,data=sub.sites,family=poisson(link=log))
gs_rs   <- glm(spec~Elev+I(Elev^2)+logpop+I(logpop^2)+yieldEVI+I(yieldEVI^2)+Transect,data=sub.sites,family=poisson(link=log))
gs_null <- glm(spec~1,data=sub.sites,family=poisson(link=log)) # nothing to be dropped

up <- stepAIC(ga_full,direction="back",trace=T)
summary(up)
plot(allEffects(up))
#(aa <- aictab(list(ga_full,ga_lu,ga_rs,ga_null,up),c("Land-use + Remote-sensing data","Land-use only","Remote-sensing only","Null model","Updated model"),second.ord = T,sort = T))
(aa <- aictab(list(ga_null,up),c("Null model","Updated model"),second.ord = T,sort = T))
#aa$rsquaredMarginal <- unlist(lapply(list(ga_full,ga_lu,ga_rs,ga_null,up),function(x)r.squaredGLMM(x)["R2m"]))
aa$LR_rsquared <- c(attr(r.squaredLR(up,ga_null),"adj.r.squared"),0)
aa <- as.data.frame(aa)
rownames(aa) <- c(up$call$formula,"Null-model")

up <- stepAIC(gs_full,direction="back",trace=T)
summary(up)
plot(allEffects(up))
#(as <-aictab(list(gs_full,gs_lu,gs_rs,gs_null,up),c("Full","LU only","RS only","Null model","Updated model"),second.ord = T,sort = T))
#as$rsquaredMarginal <- unlist(lapply(list(gs_full,gs_lu,gs_rs,gs_null,up),function(x)r.squaredGLMM(x)["R2m"]))
(as <- aictab(list(gs_null,up),c("Null model","Updated model"),second.ord = T,sort = T))
as$LR_rsquared <- c(attr(r.squaredLR(up,gs_null),"adj.r.squared"),0)
as <- as.data.frame(as)
rownames(as) <- c(up$call$formula,"Null-model")


# Make data.frame for export
aa <- aa[,c(3,4,9)];rownames(aa) <- NULL
as <- as[,c(3,4,9)];rownames(as) <- NULL
names(aa) <- c("AICc","delta_AICc","r.squareLR")
names(as) <- c("AICc","delta_AICc","r.squareLR")
rownames(aa) <- c("Land-use + Land-use intensity + Altitude","Null-model")
rownames(as) <- c("Land-use + Altitude","Null-model")

library(gridExtra)
png("OutputComparison/ModelResultsFieldwork1.png",width=8,height=6,units="in",type = "cairo",res=400)
p1 <- tableGrob(round(aa,3))
grid.arrange(p1)
dev.off()
png("OutputComparison/ModelResultsFieldwork2.png",width=8,height=6,units="in",type = "cairo",res=400)
p2 <- tableGrob(round(as,3))
grid.arrange(p2)
dev.off()

plot(resid(ga_lu)~sub.sites$Elev,col=sub.sites$Transect)
d <- data.frame(res=resid(ga_lu),elev=sub.sites$Elev,lu=sub.sites$PREDICTS.LU)
ggplot(d,aes(x=elev,y=res,group=lu,color=lu)) + geom_point() + geom_hline(yintercept=0) + facet_grid(~lu)

confint(ga_lu)
# 
inf <- influence.measures(ga_lu)
ind <- which(apply(inf$is.inf, 1, any)) # Influential points
summary(inf)
sub.sites[!inf$is.inf[,"cov.r"],]

#  Compare Cropland only of both field and PREDICTS data
sub.sites.pred <- subset(sites.pred,PREDICTS.LU%in%c("Primary Vegetation","Cropland")&PREDICTS.LUI%in%c("Minimal use","Light use"))
sub.sites.pred$PREDICTS.LU <- droplevels(sub.sites.pred$PREDICTS.LU);sub.sites.pred$PREDICTS.LUI <- droplevels(sub.sites.pred$PREDICTS.LUI)
sub.sites.pred$logpop <- log1p(sub.sites.pred$pop)# logtransform population
sub.sites <- subset(sites,PREDICTS.LU%in%c("Primary Vegetation","Cropland"))
sub.sites$PREDICTS.LU <- droplevels(sub.sites$PREDICTS.LU)
sub.sites$logpop <- log1p(sub.sites$pop) # Logtransform Population

rc <- compare_randoms(sub.sites.pred,responseVar = "logabund",fixedFactors ="CroplandSplit",otherRandoms = "SSB",fitFamily = gaussian,siteRandom = T)
pa <- glmer(paste0("logabund ~ CroplandSplit +Elev+yieldEVI+logpop+",rc$best.random),data=sub.sites.pred,family=gaussian)
rc <- compare_randoms(sub.sites.pred,responseVar = "Species_richness",fixedFactors ="CroplandSplit",otherRandoms = "SSB",fitFamily = poisson,siteRandom = T)
ps <- glmer(paste0("Species_richness ~ CroplandSplit +Elev+yieldEVI+logpop+",rc$best.random),data=sub.sites.pred,family=poisson)

#upp <- step(p)
#upp <- glmer(upp$model@call,data=sub.sites.pred,family=gaussian)
plot(allEffects(ps))
plot(Effect(focal.predictors = "CroplandSplit",pa))

fa <- glm("logabund ~ CroplandSplit",data=sub.sites,family=gaussian)
fs <- glm("spec ~ CroplandSplit",data=sub.sites,family=gaussian)
#up <- stepAIC(fa,direction="back",trace=T)
summary(fs)
plot(allEffects(fa))
plot(Effect(focal.predictors = "CroplandSplit",fa))

# Construct Data
a = (allEffects(pa)$CroplandSplit)
b = (allEffects(fa)$CroplandSplit)

pind <- which(sub.sites.pred$CroplandSplit=="Minimal use Cropland")
find <- which(sub.sites$CroplandSplit=="Minimal use Cropland")

pairs.panels(c(data.frame(resid(fs)[find]),sub.sites[find,c(55:60)]),hist.col = "grey80")

#### Include FW data and model coefficients ####
a <- data.frame(SS=sites$Transect,logabund=sites$logabund,Species_richness=sites$spec,PREDICTS.LU=sites$PREDICTS.LU)
b <- data.frame(SS=sites.pred$SS,logabund=sites.pred$logabund,Species_richness=sites.pred$Species_richness,PREDICTS.LU=sites.pred$PREDICTS.LU)
d <- rbind(a,b)

fit <- glmer(Species_richness ~ PREDICTS.LU + (1|SS),data = d,family = poisson)
plot(ranef(fit),col=c("red","red",rep("blue",nrow(ranef(fit)$SS)-2)))

#### Auto correlation in data ? ####

# Using roquefort

require(pgirmess)
library(ncf)


corr_spec <- ncf::correlog(x = sites$Lat,y=sites$Long,z = sites$spec,increment = 1,resamp = 1000,latlon = T)
corr_abund <- ncf::correlog(x = sites$Lat,y=sites$Long,z = sites$logabund,increment = 1,resamp = 1000,latlon = T)

plot.correlog(corr_spec);abline(h=0)
plot.correlog(corr_abund);abline(h=0)

ldcorr_spec <- ncf::correlog(x = sites.pred$Latitude,y=sites.pred$Longitude,z = sites.pred$logabund,increment = 1,resamp = 10,latlon = T,na.rm = T)
plot.correlog(ldcorr_spec);abline(h=0)


coords <- cbind(sites$Lat,sites$Long)
lm.morantest(gs_lu,listw = nb2listw(knn2nb(knearneigh(coords,longlat = T,k = 3)), style="W"))

moran.test(sites$spec, nb2listw(knn2nb(knearneigh(coords,longlat = T,k = 2)), style="W"))



#### Cropland residual plotting against Covariates #####
# Are my studies somehow
# Get RMSE for simple landuse model with my fielddata included
# for residuals
# plot against pop and yieldEVI
names(sites.pred)
sub <- subset(sites.pred,select = c("logabund","Species_richness","PREDICTS.LU","SS","new.mean.evi","new.yield.evi"))
sub2 <- subset(sites,select=c("logabund","spec","PREDICTS.LU","Transect","new.mean.evi","new.yield.evi"))
names(sub2) <- c("logabund","Species_richness","PREDICTS.LU","SS","new.mean.evi","new.yield.evi")
sub <- rbind(sub,sub2)

lu_ab <- glmer(logabund ~ PREDICTS.LU*new.yield.evi+ (PREDICTS.LU|SS),data=sub,family=gaussian,REML=T)
lu_sp <- glmer(Species_richness ~ PREDICTS.LU  * new.yield.evi +  (PREDICTS.LU|SS),data=sub,family=poisson(link=log),REML=T)

ag <- aggregate(sub$new.yield.evi,by=list(SS=sub$SS),FUN=function(x)mean(x,na.rm=T))
ag <- aggregate(sub$new.mean.evi,by=list(SS=sub$SS),FUN=function(x)mean(x,na.rm=T))
co <- coef(lu_ab)$SS
ag$col <- c(rep(1,60),16,16)
plot(co$PREDICTS.LUCropland~ag$x,pch=ag$col);abline(h=0)

resid(lu_sp)^2
sqrt(mean((df$model-df$measure)^2,na.rm=TRUE))

rmse <- sqrt((resid(lu_sp)^2))
plot(rmse~sub$new.yield.evi)


ag <- aggregate(sub$new.mean.evi,by=list(sub$PREDICTS.LU),FUN=function(x)mean(x,na.rm=T))

plot(fixef(lu_ab)~ag$x)
lm(fixef(lu_sp)~ag$x)



# Potential BIAS in the data and models---------------

library(influence.ME)
library(dplyr)
library(magrittr)
e_a <- exclude.influence(fit_a,obs = T)
e_alu <- exclude.influence(fit_alu,obs = T)
e_s <- exclude.influence(fit_s,obs = T)
e_slu <- exclude.influence(fit_slu,obs = T)
# No difference

# Test for sampling completeness
pred <- predictsdata
# Kick out unsuitable studies
pred <- subset(pred,Source_ID!="GP1_2012__Strauch")
pred <- subset(pred,!(Predominant_habitat%in%c("Pasture")))
pred <- subset(pred,year(pred$Sample_start_earliest)>=2000)

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

#  The Lomolino model (SSlomolino) is Asym/(1 + slope^log(xmid/area)) (Lomolino 2000, Dengler 2009).
# Parameter Asym is the asymptotic maximum number of species, 
# slope is the maximum slope of increse of richness, 
# and xmid is the area where half of the maximum richess is achieved. 
#Dengler, J. (2009) Which function describes the species-area relationship best? A review and empirical evaluation. Journal of Biogeography 36, 728–744. 
pred <- subset(pred,Study_common_taxon=="Aves")
res <- data.frame()
for( ss in unique(pred$SS)){
  print(ss)
  d <- subset(pred,SS==ss)
  d <- d[which(!is.na(d$PREDICTS.LU)),]  
  for(p in unique(d$PREDICTS.LU)){
  print(paste(ss,"-",p))
  dd <- d %>% subset(.,PREDICTS.LU==p) %>% acast(.,Site_name~Best_guess_binomial,value.var = "Measurement",fun.aggregate = sum)
  l <- tryN(fitspecaccum(round(dd,0),model = "lomolino",method = "exact",control=nls.control(warnOnly = F)))
  res <- rbind(res, data.frame(SS=ss,PREDICTS.LU=p,miss= as.numeric(unique(ifelse(is.na(l),NA,(1- ( nrow(dd) / predict(l,coef(l)[1]) )) * 100 ) )),N=nrow(dd)) )
  }
  #rm(d,dd,l,p)
}
#res %>% group_by(PREDICTS.LU) %>% summarise(AvM = mean(miss,na.rm=T))
d <- species %>% melt(.) %>% merge(.,sites,by.x="Var1",by.y="SiteName") 
res2 <- data.frame()
for(p in unique(d$PREDICTS.LU)){
  print(p)
  dd <- d %>% subset(.,PREDICTS.LU==p) %>% acast(.,Var1~Var2,value.var = "value",fun.aggregate = sum)
  l <- tryN( fitspecaccum(dd,model = "lomolino",method = "exact"))
  res2 <- rbind(res2, data.frame(SS=ss,PREDICTS.LU=p,miss= as.numeric(unique(ifelse(is.na(l),NA,(1- ( nrow(dd) / predict(l,coef(l)[1]) )) * 100 ) )),N=nrow(dd) ) )
}

#res2 %>% group_by(PREDICTS.LU) %>% summarise(AvM = mean(miss,na.rm=T))
rm(l,d,dd,p,ss)
res$Data <- "Broad-scale bird studies \n (Mean accross studies)"
res2$Data <- "Independent"
r <- rbind(res,res2)
r$miss <- 1-(r$miss * 0.01) # covnert back to percentages
r <- remNa(r)
r <- ddply(r,.(PREDICTS.LU,Data),summarise,miss=mean(miss,na.rm = T),N=paste0("N=",sum(N)) )
r$PREDICTS.LU <- factor(r$PREDICTS.LU,levels = ord)
rect_colours<-c("#00AE00","#94BD5E", "#006B6B", "#E6E64C","#808080")
rect_lab <- c("PV","SV","PL","CL","UR")

ggthemr("fresh",layout = "clean",type = "outer",spacing = 0.5)
g <- ggplot(r,aes(x=PREDICTS.LU,y=miss,fill=PREDICTS.LU,group=factor(Data)) )
g <- g + geom_bar(stat = "identity",position="dodge")
g <- g + geom_hline(yintercept=1,color="red", lty=2,size=1.2,alpha=.6)
g <- g + scale_x_discrete(labels=rect_lab)
g <- g + facet_grid(~Data,scales = "fixed")
g <- g + geom_text(aes(label=N,size=15),color="black",position=position_dodge(width=0.9),hjust = 0.5, vjust = 2,show_guide=F)
g <- g + scale_y_continuous(labels = percent_format(),breaks=pretty_breaks(),expand=c(0,0))
g <- g + xlab("") + ylab("Percentual distance to Assymptote (100%)") # No labels  
g <- g + scale_fill_manual(values=rect_colours,guide = "none") 
g <- g + theme(axis.title.x=element_blank(),axis.ticks.x=element_blank()) # Remove legend and x-axis stuff
g <- g + ggtitle("Sampling completeness")
g

#FOR THESIS
ggsave("../000Thesis_writeup/WriteUpLatex/gfx/SIFigure3.pdf",plot=g+theme(legend.position="none"),scale=1.1,dpi=400)

#--- Export table for lookup-----------
pred <- predictsdata
# Kick out unsuitable studies
pred <- subset(pred,Source_ID!="GP1_2012__Strauch")
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


ddply(pred,.(Source_ID),summarise,
      N=n(),
      Per=(n()/nrow(pred))*100 )

p <- subset(pred,select = c("SS","Site_name","Longitude","Latitude","LandUse"))
p <- unique(p)

write.csv2(p,"~/o.csv",row.names=F)
