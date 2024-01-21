###############################################################################
#       Modeling and Mapping Adaptability of Eucalyptus Clones                #
#                           version 2.0  - English                            #
#                                                                             #
# Elaborado por:Leonardo Oliveira Silva da Costa                              #
#                                                                             #
###############################################################################


#Packages and functions:
rm(list=ls())

library(raster)
library(terra)
library(ggplot2)
library(gridExtra)
library(dplyr)
source("Functions_GIS_PLS.R")

############################################################################
####################                                  ######################
##################        Importing and Configuring     ####################
##################         Inventory Table              ####################
####################                                  ######################
############################################################################

inv<- read.csv("mai_set.csv", header=T,sep=';');inv$id<-1:nrow(inv)
str(inv)

#Date/Year
inv$planting_date<- as.Date(inv$planting_date, "%d/%m/%Y")
inv$year<- format(inv$planting_date,"%Y")
inv$mon<- format(inv$planting_date,"%m")

#Semester
p<-paste0("0",1:6)
s<-c(paste0("0",7:9),10:12)

inv <- inv %>% mutate(sem = case_when(mon %in% p ~ 1,mon %in% s ~2))

#Year.Semester
inv$year.sem<- paste(inv$year, inv$sem, sep=".")

############################################################################
####################  Associating Each Pixel Value with Data  ##############
############################################################################

# Downloading 5-minute raster for a 10 km² grid template
getWC(var="elev",res="r5min") #Molde de 5 minutos

#Importando raster molde
raster<- raster("Covariates/r5min/wc2.1_5m_elev.tif")
plot(raster) #Conferindo

# Transforming inventory table into points
coordinates(inv)<- ~ lon + lat 
plot(inv,add=T)

# Obtaining points for each grid
pixel<-grid_point(raster = raster, point = inv)
pixel<- merge(data.frame(inv),pixel) # Merging inventory data and pixel data


# Averaging by quadrant and semester

(data<-pixel %>% group_by(gen, ID_Pixel,reg,year.sem) %>% summarise( mai = mean(mai),
                                                                               lat = mean(lat),
                                                                               lon = mean(lon),
                                                                               year = mean(as.numeric(year)),
                                                                               sem = mean(sem),
                                                                               datai= min(planting_date)))

nrow(data)

write.table(data,"pixel_5m_toyset.csv",sep=';')
############################################################################
##################                                     #####################
####################     Environmental Typing        #######################
#######################                          ###########################
############################################################################

#Download rasters

getWC(var=c("bioc","elev"), res="r2.5min") #WorldClim (Biovars e Elevation)
getSG(attributes = c("bdod","clay","sand","cec"),resample = T,raster =raster("./Covariates/r2.5min/wc2.1_2.5m_elev.tif") )
getLonLat(raster = raster("./Covariates/r2.5min/wc2.1_2.5m_elev.tif"))

# Cutting the rasters for the area of interest
e<-stack_rasters(path = "./Covariates/r2.5min/",pattern = ".tif")
plot(e)

br<-shapefile("./Shapefiles/BrasilWGS.shp") # Importing the template shapefile
e<-map_cut(e,br) # Cutting all the rasters
plot(e)

# Saving the stack
writeRaster(e,"./Covariates/r2.5m_Covariates_Cut.tiff")

#Data
coordinates(data) <- ~ lon + lat
names(e)

# The rastertype function extracts the covariates from the rasters for each point
df.cov <- rastertype(stack = e[[3:26]], points = data) # Excluding lat and lon
str(df.cov)

write.table(df.cov,"Covariates_5m_toyset.csv",sep=';')

############################################################################
##################                                     #####################
####################           GGEPLS              #########################
#######################                          ###########################
############################################################################

data_cov<-read.csv("Covariates_5m_toyset.csv",sep=';',header=T)


str(data_cov)
str(data_cov[c(6:7,12:35)])

config<-list(env="reg", genotype= "gen", pred = c(6:7,12:35), resp=5)

fit<-fit_ggepls(dados = data_cov,
                genotype = "gen",
                env  = "reg",resp =5, pred = c(6:7,12:35))
#or

fit<-fit_ggepls(dados = data_cov,config = config)

fit$metrics
fit$values
fit$coef
fit$std_coef
fit$circle_coord
fit$response
fit$predictors 
fit$genotypes


############################################################################
##################                                     #####################
####################           Validation          #########################
#######################                          ###########################
############################################################################

loo<-validation_ggepls(data_cov,config = config)

loo$metrics
loo$values

############################################################################
##################                                     #####################
#############            Correlation Circles             ###################
#######################                          ###########################
############################################################################

fit$coef
ca<-covariates_analysis(ggepls = fit, rank=5,resp_name = "MAI",pred_names=c("LATI","LONG","BDOD","CEC","CLAY","SAND",paste0("BIO",c(1,10:19,2:9)),"ELEV"))

ca$ranking
ca$cicles[[4]] #1 to genotype number
ca$frequency

############################################################################
##################                                     #####################
####################            Maps                ########################
#######################                          ###########################
############################################################################


# Standardizing the names of the stack with the coefficients
fit$coef
names(e)
names(e)[1:2]<-c("lat","lon") 

coordinates(data) <- ~ lon + lat
buf<-shapefile("./Shapefiles/Alvo_100km.shp")

maps<-ggepls_map(ggepls = fit,stack = e);plot(maps) # Plotting for the entire region of the rasters
maps_cut<-ggepls_map(ggepls = fit,stack = e,cut=T,shapefile = buf);plot(maps_cut)  # Plotting after cutting

# It's also possible to cut later without the cut option, using the separate function
maps_cut2<-map_cut(maps,buf)
plot(maps_cut2)

############################################################################
##################                                     #####################
################          Which Won Where               ####################
#######################                          ###########################
############################################################################

# Which won where without restrictions
w<-which_won_where(maps)

plot_recomendation(w)

# Which won where:
# Removing clone CLZ003
# By region of the shapefile "buf"

# Create a matrix of clone presence by region
(X<-get_X(data_cov,genotypes<-"gen", region<-'reg'))

# Clone recommendation among those planted in each region
w<-which_won_where(maps,rm_gen = c("CLZ003"),by_region = T,shapefile = buf,regions = "unf",X = X)

w$general_recomendation
w$region_recomendation

#Plot 
(g<-plot_recomendation(w))
plot(w$raster)

