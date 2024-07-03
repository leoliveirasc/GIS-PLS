###############################################################################
#       Modelagem e Mapeamento de adaptabilidade de clones de eucalipto       #
#                                  Versao 2.0  - PT BR                       #
#                                                                              #
# Elaborado por:Leonardo Oliveira Silva da Costa                              #
#                                                                             #
###############################################################################

#Pacotes e funcoes:
rm(list=ls())

library(raster)
library(terra)
library(ggplot2)
library(gridExtra)
library(dplyr)
source("./Functions_GISPLS/Functions_GIS_PLS.R")

############################################################################
####################                                  ######################
##################        Importando e configurando     ####################
##################        tabela de inventário          ####################
####################                                  ######################
############################################################################

inv<- read.csv("mai_set.csv", header=T,sep=';');inv$id<-1:nrow(inv)
str(inv)

#DATA/ANO
inv$planting_date<- as.Date(inv$planting_date, "%d/%m/%Y")
inv$year<- format(inv$planting_date,"%Y")
inv$mon<- format(inv$planting_date,"%m")

#SEMESTRE
p<-paste0("0",1:6)
s<-c(paste0("0",7:9),10:12)

inv <- inv %>% mutate(sem = case_when(mon %in% p ~ 1,mon %in% s ~2))

#ANO.SEMESTRE
inv$year.sem<- paste(inv$year, inv$sem, sep=".")

############################################################################
#################### Associando cada valor de pixel aos dados ##############
############################################################################

#Download de raster 5min para molde de grid 10km²
getWC(var="elev",res="r5min") #Molde de 5 minutos

#Importando raster molde
raster<- raster("Covariates/r5min/wc2.1_5m_elev.tif")
plot(raster) #Conferindo

#transformando tabela de inventário em pontos
coordinates(inv)<- ~ lon + lat 
plot(inv,add=T)

#Obtendo os pontos de cada grid
pixel<-gridPoint(raster = raster, point = inv)
pixel<- merge(data.frame(inv),pixel) #Juntando os dados de inventário e os dados 

#Médias a partir do quadrante e semestre

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
####################        Tipagem ambiental        #######################
#######################                          ###########################
############################################################################

#Download rasters

getWC(var=c("bioc","elev"), res="r2.5min") #WorldClim (Biovars e Elevacao)
getSG(attributes = c("bdod","clay","sand","cec"),resample = T,raster =raster("./Covariates/r2.5min/wc2.1_2.5m_elev.tif") )
getLonLat(raster = raster("./Covariates/r2.5min/wc2.1_2.5m_elev.tif"))

#Cortando os rasters para a area de interesse
e<-stackRasters(path = "./Covariates/r2.5min/",pattern = ".tif")
plot(e)

br<-shapefile("./Shapefiles/BrasilWGS.shp") #importando o shapefile de molde
e<-mapCut(e,br) #cortando todos os rasters
plot(e)

#Salvando o stack
writeRaster(e,"./Covariates/r2.5m_Covariates_Cut.tiff")

#Dados

coordinates(data) <- ~ lon + lat
names(e)

#funcao rastertype extrai as covariaveis dos rasters para cada ponto
df.cov<-rastertype(stack = e[[3:26]], points = data) #Excluindo lat e lon 

str(df.cov)
write.table(df.cov,"Covariates_5m_toyset.csv",sep=';')

############################################################################
##################                                     #####################
####################           Ajuste              #########################
#######################                          ###########################
############################################################################

#data<-read.csv("pixel_5m_toyset.csv",sep=';',header=T)
data_cov<-read.csv("Covariates_5m_toyset.csv",sep=';',header=T)

#data_cov<-cbind(data.frame(data),df.cov)
str(data_cov)
str(data_cov[c(6:7,12:35)])

config<-list(env="reg", genotype= "gen", pred = c(6:7,12:35), resp=5)

fit<-fitGISPLS(data = data_cov,
                genotype = "gen",
                env  = "reg",resp =5, pred = c(6:7,12:35))
#OU

fit<-fitGISPLS(data = data_cov,config = config)

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
####################           Validacao           #########################
#######################                          ###########################
############################################################################

loo<-validateGISPLS(data_cov,config = config)

loo$metrics
loo$values

############################################################################
##################                                     #####################
#############          Circulos de correlacao            ###################
#######################                          ###########################
############################################################################

fit$coef
ca<-covariatesAnalysis(gispls = fit, rank=5,resp_name = "MAI",pred_names=c("LATI","LONG","BDOD","CEC","CLAY","SAND",paste0("BIO",c(1,10:19,2:9)),"ELEV"))

ca$ranking
ca$cicles[[4]] #1 a n de genótipos
ca$frequency

############################################################################
##################                                     #####################
####################           Mapas                ########################
#######################                          ###########################
############################################################################


#Padronizando os nome do stack com os dos coeficientes
fit$coef
names(e)
names(e)[1:2]<-c("lat","lon") 

coordinates(data) <- ~ lon + lat
buf<-shapefile("./Shapefiles/Alvo_100km.shp")

maps<-mapGISPLS(gispls = fit,stack = e);plot(maps) #Plotagem para toda a regiao dos rasters
maps_cut<-mapGISPLS(gispls = fit,stack = e,cut=T,shapefile = buf);plot(maps_cut) #Plotagem com corte

#tbm eh possivel cortar depois sem a opcao cut, utilizando a funcao separada
maps_cut2<-mapCut(maps,buf)
plot(maps_cut2)

############################################################################
##################                                     #####################
################          Quem vence onde             ######################
#######################                          ###########################
############################################################################

##Quem vence onde sem restricoes
w<-whichWonWhere(maps)
recommendationPlot(w)

#wich won where:
# retirando clone CLZ003
# Por regiao do shapefile br

#Criar matriz de presenca dos clones por regiao
(X<-getX(data_cov,genotypes<-"gen", region<-'reg'))

#Recomendacao de clones dentre aqueles plantados em cada regiao
w<-whichWonWhere(maps,rm_gen = c("CLZ003"),by_region = T,shapefile = buf,regions = "unf",X = X)

w$general_recomendation
w$region_recomendation

#Plot 
(g<-recommendationPlot(w))
plot(w$raster)
