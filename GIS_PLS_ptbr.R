###############################################################################
#       Modelagem e Mapeamento de adaptabilidade de clones de eucalipto       #
#                                  VersÃ£o 2.0                                #
#                                                                             #
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
source("Functions_GIS_PLS.R")
source("cov_group_selection.R")

############################################################################
####################                                  ######################
##################        Importando e configurando     ####################
##################        tabela de inventário          ####################
####################                                  ######################
############################################################################

inv<- read.csv("dados_inventario_filtrado.csv", header=T,sep=';');inv$id<-1:nrow(inv)
str(inv)

#DATA/ANO
inv$data_plantio<- as.Date(inv$data_plantio, "%d/%m/%Y")
inv$data_medicao<- as.Date(inv$data_medicao, "%d/%m/%Y")
inv$ANOi<- format(inv$data_plantio,"%Y")
inv$ANOf<- format(inv$data_medicao,"%Y")
inv$MESi<- format(inv$data_plantio,"%m")
inv$MESf<- format(inv$data_medicao,"%m")

#SEMESTRE
p<-paste0("0",1:6)
s<-c(paste0("0",7:9),10:12)

inv <- inv %>% mutate(SEMi = case_when(MESi %in% p ~ 1,MESi %in% s ~2),
                      SEMf = case_when(MESf %in% p ~ 1,MESf %in% s ~2))

#ANO.SEMESTRE
inv$ANO.SEM<- paste(inv$ANOi, inv$SEMi, sep=".")

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
pixel<-grid_point(raster = raster, point = inv)
pixel<- merge(data.frame(inv),pixel) #Juntando os dados de inventário e os dados 

#Médias a partir do quadrante e semestre

(data<-pixel %>% group_by(MatGen_oficial, ID_Pixel,unf,ANO.SEM) %>% summarise( ima7cc = mean(ima7cc),
                                                                              lat = mean(lat),
                                                                              lon = mean(lon),
                                                                              ANOi = mean(as.numeric(ANOi)),
                                                                              SEMi = mean(SEMi),
                                                                              idade = mean(idade_inv),
                                                                              dataf= max(data_medicao),
                                                                              datai= min(data_plantio),
                                                                              alt = mean(altitude_raw)))
nrow(data)

write.table(data,"pixel_5m.csv",sep=';')
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
e<-stack_rasters(path = "./Covariates/r2.5min/",pattern = ".tif")
plot(e)

br<-shapefile("./Shapefiles/BrasilWGS.shp") #importando o shapefile de molde
e<-map_cut(e,br) #cortando todos os rasters
plot(e)

#Salvando o stack
writeRaster(e,"./Covariates/r2.5m_Covariates_Cut.tiff")

#Dados

coordinates(data) <- ~ lon + lat
names(e)

#funcao rastertype extrai as covariaveis dos rasters para cada ponto
df.cov<-rastertype(stack = e[[3:26]], points = data) #Excluindo lat e lon 

str(df.cov)
write.table(df.cov,"Covariates_5mauto.csv",sep=';')

############################################################################
##################                                     #####################
####################           Ajuste              #########################
#######################                          ###########################
############################################################################

data<-read.csv("pixel_5m.csv",sep=';',header=T)
df.cov<-read.csv("Covariates_5mauto.csv",sep=';',header=T)

data_cov<-cbind(data.frame(data),df.cov)
str(data_cov)
str(data_cov[c(6:7,28:51)])

config<-list(env="unf", genotype= "MatGen_oficial", pred = c(6:7,28:51), resp=5)

fit<-fit_ggepls(dados = data_cov,
                genotype = "MatGen_oficial",
                env  = "unf",resp =5, pred = c(6:7,28:51))
#OU

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
####################           Validacao           #########################
#######################                          ###########################
############################################################################

loo<-validation_ggepls(data_cov,config = config)

loo$metrics
loo$values

############################################################################
##################                                     #####################
#############          Circulos de correlacao            ###################
#######################                          ###########################
############################################################################


ca<-covariates_analysis(ggepls = fit, rank=5,resp_name = "MAI",pred_names=c("LATI","LONG","BDOD","CEC","CLAY","SAND",paste0("BIO",1:19),"ELEV"))

ca$ranking
ca$cicles[[3]] #1 a n de genótipos
ca$frequency

############################################################################
##################                                     #####################
####################           Mapas                ########################
#######################                          ###########################
############################################################################

fit$coef
names(e)[1:2]<-c("lat","lon") # Padronizando os nome do stack com os dos coeficientes

coordinates(data) <- ~ lon + lat
buf<-shapefile("./Shapefiles/Alvo_100km.shp")

maps<-ggepls_map(ggepls = fit,stack = e);plot(maps) #Plotagem para toda a regiao dos rasters
maps_cut<-ggepls_map(ggepls = fit,stack = e,cut=T,shapefile = buf);plot(maps_cut) #Plotagem com corte

#tbm eh possivel cortar depois sem a opcao cut, utilizando a funcao separada
maps_cut2<-map_cut(maps,buf)
plot(maps_cut2)

############################################################################
##################                                     #####################
################          Quem vence onde             ######################
#######################                          ###########################
############################################################################

##wich won where sem restricoes
w<-which_won_where(maps)

plot_recomendation(w)

#wich won where:
# retirando clone CLZ003
# Por regiao do shapefile br

#Criar matriz de presenca dos clones por regiao
(X<-get_X(data_cov,genotypes<-"MatGen_oficial", region<-'unf'))

#Recomendacao de clones dentre aqueles plantados em cada regiao
w<-which_won_where(maps,rm_gen = c("CLZ003"),by_region = T,shapefile = buf,regions = "unf",X = X)

w$general_recomendation
w$region_recomendation

#Plot com a legenda correta
(g<-plot_recomendation(w))
plot(w$raster)
w

