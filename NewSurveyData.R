###########Data prep############
library(readxl)
library(dplyr)
library(marfunky)
detachPackage("MODIS");detachPackage("raster")
dataSur <- read_excel("../FieldData_Surveys/BirdButtCollate.xlsx",1)
source("Validation.R")
valid_column_names <- make.names(names=names(sites), unique=TRUE, allow_ = TRUE)
names(sites) <- valid_column_names
sitesD <- sites %>% filter(Authority=="DickensO",Transect=="Taita")

# Select and aggregate per species and plot
dataSur <- as.data.frame(dataSur)
names(dataSur) <- make.names(names(dataSur))
sm <- acast(dataSur,PLOT.NO~BINOMIAL,fun.aggregate = sum,value.var="No..Indiv")
sms <- data.frame(PLOT.NO = row.names(sm),
                  RSabund = apply(sm,1,function(x) sum(x,na.rm = T)),
                  RSspec = rowSums(sm > 0,na.rm = T) )

sms <- dataSur %>% dplyr::select(PLOT.NO,Latitude,Longitude) %>%
  mutate(PLOT.NO = as.factor(PLOT.NO)) %>% unique() %>% 
  right_join(.,sms) %>% mutate(RSlogabund = log1p(RSabund))

sms$Long_r <- round(sms$Longitude,6)
sms$Lat_r <- round(sms$Latitude,6)

# Prep for merging
sitesD$Long_r <- round(sitesD$Longitude,6)
sitesD$Lat_r <- round(sitesD$Latitude,6)
sitesD$PLOT.NO <- numbers_from_string(sitesD$SiteName)
a <- join(sitesD,sms,by=c("PLOT.NO" = "PLOT.NO"))

qplot(a$spec,a$RSspec,shape=a$PREDICTS.LU,size=5)

setdiff2(sitesD$Long_r,dataSur$Long_r) # 4 resurveyed Sites not in the original
setdiff2(sitesD$Lat_r,dataSur$Lat_r) # 4 resurveyed Sites not in the original
match(sitesD$Long_r,dataSur$Long_r)

numbers_from_string(sitesD$SiteName)
j <- right_join(dataSur,sitesD)
