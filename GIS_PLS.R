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
library(plsdepot)
library(ggplot2)
library(gridExtra)
library(dplyr)
source("Functions_GIS_PLS.R")

############################################################################
####################                                  ######################
##################        Importando e configurando     ####################
##################        tabela de inventário          ####################
####################                                  ######################
############################################################################

inv<- read.csv("./Tabelas/dados_inventario_filtrado.csv", header=T,sep=';');inv$id<-1:nrow(inv)
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

#Importando raster molde
raster<- raster("./Rasters_5m/wc2.1_5m_elev.tiff")

#transformando tabela de inventário em pontos
coordinates(inv)<- ~ lon + lat 

#Obtendo os pontos de cada grid
pixel<-grid_point(raster = raster, point = inv, point_id = "id")
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

#Funçao para juntar todos os rasters em um só

e<-stack_rasters(path = "./Rasters_2.5m/",pattern = ".tif" )

#Função que extrai os pontos
#Colocar um arquivo intermediário para limpar a memoria

coordinates(data) <- ~ lon + lat
df.cov<-rastertype(stack = e[[1:25]], points = data)
str(df.cov)
write.table(df.cov,"Covariates_5m.csv",sep=';')

############################################################################
##################                                     #####################
####################           Ajuste              #########################
#######################                          ###########################
############################################################################

data<-read.csv("pixel_5m.csv",sep=';',header=T)
df.cov<-read.csv("Covariates_5m.csv",sep=';',header=T)
data_cov<-cbind(data.frame(data),df.cov)
str(data_cov)
data_cov<-data_cov[-18] #Tirar Silte


config<-list(env= "unf", genotype= "MatGen_oficial", pred = c(6:7,14:37), resp=5)

fit<-fit_ggepls(dados = data_cov,
                genotype = "MatGen_oficial",
                env  = "unf",resp =5, pred = c(6:7,14:37))
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

names(data_cov[config$pred])

############################################################################
##################                                     #####################
####################           Circulos            #########################
#######################                          ###########################
############################################################################


ca<-covariates_analysis(ggepls = fit, rank=5,resp_name = "MAI",pred_names=c("LATI","LONG","BDOD","CEC","CLAY","SAND",paste0("BIO",1:19),"ELEV"))

ca$ranking
ca$cicles[[1]] #1 a n de genótipos

############################################################################
##################                                     #####################
####################           Mapas                ########################
#######################                          ###########################
############################################################################

e<-stack_rasters(path='./Rasters_2.5m/',pattern = '.tif')
names(e)
names(e)[25:26]<-c("lat","lon") # Padronizando os nome do stack com os dos coeficientes
br<-shapefile("./Shapefiles/Alvo_100km.shp") #Shapefile molde para corte

maps<-ggepls_map(ggepls = fit,stack = e,cut=T,shapefile = br)

plot(maps)

#tbm eh possivel cortar depois sem a opcao cut, utilizando a funcao separada
maps_cut<-map_cut(maps,br)
plot(maps_cut)

############################################################################
##################                                     #####################
################          Quem vence onde             ######################
#######################                          ###########################
############################################################################

(X<-get_X(data_cov,genotypes<-"MatGen_oficial", region<-'unf'))

w<-which_won_where(maps,rm_gen = "CLZ003",by_region = T,shapefile = br,regions = "unf",X = X)  
  
plot(w,col=rainbow(6))  

