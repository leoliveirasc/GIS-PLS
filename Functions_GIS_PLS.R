
source("plsreg1_mod.R")

##########################################################
######## Grid para obter as medias por quadrante #########
##########################################################

grid_point<-function(raster,point){
  #point<-inv; raster<- raster; point_id<-"id"
  
  values(raster)<-1:length(values(raster))
  
  #Barra de progresso
  total <- nrow(point@data)
  print("Obtendo os pontos de cada grid")
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  for (i in 1:nrow(point@data)){
    
    if(i==1) vars<- vector()
    
    indv<- point[i,]
    pp<- raster::extract(raster,indv)
    vars<- c(vars, pp)
    
    #Atualizando Barra de Progresso
    setTxtProgressBar(pb, i)
    
  }
  
  
  mt<-matrix(vars, nrow=nrow(point), ncol=1)
  df.pixel<- data.frame(mt); names(df.pixel)<-"ID_Pixel"
  df.pixel$id<- point$id
  
  return(df.pixel)
  
}
##########################################################
##################  Download WorldClim  ##################
##########################################################

getWC<-function(var=NULL,res=NULL){
  #var='elev';res="r5min"
  links<-read.csv("word.links.csv",h=T,sep=';')
  
  if("./Covariates" %in% list.dirs() == F) dir.create("Covariates")
  if(paste0("./Covariates/",res) %in% list.dirs() == F) dir.create(paste0("./Covariates/",res))
  
  file<-links[(links$variable%in%var) &(links$res == res),]$link

  options(timeout=6000)
  
  for(i in 1:length(file)){
  download.file(destfile = "./Covariates/getWC.zip",url = file[i])
  unzip("./Covariates/getWC.zip",exdir = paste0("./Covariates/",res,"/"))
  }
  options(timeout=60)
  file.remove("./Covariates/getWC.zip")
  
}
##########################################################
##################  Download SoilGrids  ##################
##########################################################

library('XML')
library(gdalUtilities)
library(raster)

#Profundidades
#0-5cm_mean
#5-15cm_mean
#15-30cm_mean
#30-60cm_mean
#60-100cm_mean
#100-200cm_mean


#Propriedades
#bdod
#cec
#cfvo
#clay
#nitrogen
#ocd
#ocs
#phh2o
#sand
#silt
#soc - soil organic carbon 
#wrb


#vetores com os atributos e profundidades

getSG<-function(attributes = NULL, layers ="30-60cm_mean", resample = F, raster=NULL,tr=c(4500,4500)){
  
  #attributes = "bdod";layers ="30-60cm_mean";tr=c(4500,4500);resample = T;raster=raster("./Covariates/r2.5min/wc2.1_2.5m_elev.tif")
    
    if("./Covariates" %in% list.dirs() == F) dir.create("Covariates")
    if("./Covariates/getSG" %in% list.dirs() == F) dir.create("./Covariates/getSG")
  
  sg_url="/vsicurl/https://files.isric.org/soilgrids/latest/data/"
  
  if(resample ==F){
    for (i in 1:length(attributes)) {
      for (j in 1:length(layers)) {
        cat("\nDOWNLOADING",paste0(attributes[i],layers[j],'\n'))
        
        datos = paste0(attributes[i], "/", attributes[i], "_", layers[j], ".vrt") #attr + layer
        lfile = paste0("./Covariates/getSG/sg.", attributes[i], layers[j],".tiff") #output path
        
        gdal_translate(paste0(sg_url, datos), lfile,tr =tr)
        
        
      }
    }
  }else{
    for (i in 1:length(attributes)) {
      for (j in 1:length(layers)) {
        cat("\n\nDOWNLOADING:",paste0(attributes[i],layers[j],'\n'))
        
        datos = paste0(attributes[i], "/", attributes[i], "_", layers[j], ".vrt") #attr + layer
        lfile = paste0("./Covariates/getSG/sg.", attributes[i], layers[j],".tiff") #output path
        
        gdal_translate(paste0(sg_url, datos), lfile,tr = tr)
        
        #Changing the projection and ressampling
        cat("\nRessampling...\n")
        gs<-rast(lfile)
        gs<-project(gs,projection(raster))
  
        gs.res<-resample(gs,terra::rast(raster))
        res<- paste0("r",res(gs.res)[1]*60,"min")
        
        if(paste0("./Covariates/",res) %in% list.dirs() == F) dir.create(paste0("./Covariates/",res))
        writeRaster(gs.res,paste0("./Covariates/",res,"/sg.", attributes[i], layers[j],".tiff"), overwrite=TRUE)
        ?writeRaster
        cat("DONE!")
 

      }
    }
  }
}


#########################################################
############# Construindo rasters Lon e Lat #############
#########################################################

getLonLat<-function(raster){
  
  #importando raster molde

  #Criando rasters auxiliares para lat e long
  lat<-lon<-raster
  
  #Atribuindo os valores para cada um dos rasters
  values(lat)<- yFromCell(raster,cell = 1:ncell(raster))
  values(lon)<- xFromCell(raster,cell = 1:ncell(raster))

  #Cortando para o raster molde
  lat<-mask(lat,raster)
  lon<-mask(lon,raster)
  
  if("./Covariates" %in% list.dirs() == F) dir.create("Covariates")
  res<- paste0("r",res(raster)[1]*60,"min")
  
  if(paste0("./Covariates/",res) %in% list.dirs() == F) dir.create(paste0("./Covariates/",res))
  
  writeRaster(lat,paste0("./Covariates/",res,"/latitude.tiff"), overwrite=TRUE)
  writeRaster(lon,paste0("./Covariates/",res,"/longitude.tiff"), overwrite=TRUE)
  
  return(stack(lon,lat))
  
}

##########################################################
######### Stack de rasters em uma mesma pasta ############
##########################################################

stack_rasters<-function(path,pattern){
  #path = "./Rasters_2.5m/";pattern = ".tif"
  
  (rasters<- c(list.files(path = path,pattern = pattern)))
  
  for (i in 1:length(rasters)){
    a.sub <- raster(paste0(path, rasters[i]))
    if(i == 1) { e <- a.sub } else { e <- stack(e,a.sub)}}
  
  return(e)
  
}

##########################################################
####### Extração de dados em rasters por ponto ###########
##########################################################

rastertype<- function(stack,points,point_id=NULL){
  #stack = e[[1:24]]; points = data_coord
  
  vars<-vector()
  cat("EXTRACAO DE DADOS AMBIENTAIS A PARTIR DE RASTERS\n")
  
  for (j in 1:length(names(stack))){   
    
    a<- stack[[j]]
    
    #Barra de progresso
    total <- nrow(points@data)
    cat("\nRASTER:",names(a),"\t",j,"/",length(names(stack)),"\n")
    
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    
    for ( i in 1:nrow(points@data)){
      
      indv<- points[i,]
      pp<- raster::extract(a,indv) 
      vars<- c(vars, pp)
      
      #Atualizando Barra de Progresso
      setTxtProgressBar(pb, i)
    }       
    cat("\n")
  }
  
  mt<-matrix(vars, nrow=nrow(points), ncol=length(names(stack)))
  edacli<- data.frame(mt); names(edacli)<- names(stack)
  
  if(is.null(point_id) ==T) point_id<-1:nrow(points@data)
  edacli$id<- point_id
  
  edacli<-cbind(data.frame(points),edacli)
  
  
  return(edacli)
}

##########################################################
############### Ajuste do modelo PLS  ####################
##########################################################

fit_ggepls<-function(dados,config=NULL,genotype=NULL,
                     env=NULL,resp=NULL, pred=NULL,c=2,
                     selected_cov= NULL){
  
  #config=config; dados = ;selected_cov= gcov
  
  if(is.null(config)== FALSE){
    genotype<- config$genotype
    env<- config$env
    resp<-config$resp
    pred<- config$pred
  }
  
  na<-nrow(dados)- nrow(na.exclude(dados))
  if(na!=0){
    warning(paste0("\nDroping ",na," rows with NA") )
    dados<-na.exclude(dados)
  }
  
  names(dados)[match(genotype,names(dados))]<-"genotype"
  names(dados)[match(env,names(dados))]<-"env"
  gen<-as.character(unique(dados$genotype))
  dados$id<-1:nrow(dados)
  
  for (j in 1:length(gen)){
    
    if(j==1){
      r2.leave<-vector()
      ima.pred<-vector()
      predict<-vector()
      observ<-vector()
      geno<-vector()
      env<-vector()
      gg.ima<-list()
      lat<-vector()
      point<-vector()
      
    }
    
    cat("\n\nINICIANDO AJUSTE GGEPLS:",gen[j],'\n')
    
    df.gen<-dados[dados$gen == gen[j],] 
    predictors<-df.gen[pred]
    response <- df.gen[resp]
    
    if(is.null(selected_cov)==T){
      
      total <- nrow(df.gen)
      # Barra de progresso
      pb <- txtProgressBar(min = 0, max = total, style = 3)
      
      #PLS
      pls <-plsreg1(predictors = predictors,response = response,comps = c,crosval = T)
      
      #Latent variables
      coo<-as.data.frame(pls$cor.xyt[,c('t1','t2')])
      coo$var<-rownames(coo)
      coo$gen <-gen[j]
      circle <- coo
      
    }else{
      
      gr<-selected_cov$selected[selected_cov$selected$gen == gen[j],]$groups
      sel.covs<-selected_cov$groups[selected_cov$groups$groups %in% strsplit(gr,split=';')[[1]],]$covs
      
      total <- nrow(df.gen)
      # Barra de progresso
      pb <- txtProgressBar(min = 0, max = total, style = 3)
      
      #PLS
      pls <-plsreg1(predictors = predictors[sel.covs],response = response,comps = c,crosval = T)
      
      pls$std.coefs
      #Latent variables
      coo<-as.data.frame(pls$cor.xyt[,c('t1','t2')])
      coo$var<-rownames(coo)
      coo$gen <-gen[j]
      circle <- coo
      
    }
    
    geno<-c(geno, as.character(df.gen$genotype))
    env<-c(env,df.gen$env)
    point<-c(point,df.gen$id)
    
    #Atualizando Barra de Progresso
    setTxtProgressBar(pb, total)
    
    pls.coefs<-rep(0,length(pred)+1)
    names(pls.coefs)<-c("Intercept",names(dados[pred]))
    pls.coefs_std<-pls.coefs[-1]
    
    pls.coefs[names(pls$reg.coefs)]<-pls$reg.coefs
    pls.coefs_std[names(pls$std.coefs)]<-pls$std.coefs
    
    value<-as.numeric(data.matrix(cbind(1,predictors))%*%(pls.coefs))
    predict<- c(predict,value)
    observ<- c(observ,response[,1])
    
    if(j==1) {
      beta<-pls.coefs
      beta_std<-pls.coefs_std
      CIRCLE<-circle
    }else
    {
      beta<-cbind(beta,pls.coefs)
      beta_std<- cbind(beta_std,pls.coefs_std)
      CIRCLE<- rbind(CIRCLE,circle)
      row.names(CIRCLE)<-1:nrow(CIRCLE)
    }
    
    
    if(j == length(gen)) df.loo<- data.frame(Genotype= geno, Env = env, Predicted = predict, Observed = observ)
  }
  
  require(dplyr)
  beta<- as.data.frame(beta)
  beta_std<- as.data.frame(beta_std)
  names(beta)<-names(beta_std)<-gen
  df.metrics<-df.loo %>% group_by(Genotype) %>% summarise(R2 = round(cor(Predicted, Observed) ^ 2,4),
                                                          RMSE = sqrt(mean((Predicted - Observed)^2)))
  
  
  return(list(values = df.loo,metrics = df.metrics,
              coef= beta, std_coef = beta_std,
              circle_coord= CIRCLE, is.ggepls =TRUE,
              response= names(response), predictors = names(predictors), genotypes =gen))
}


#############################################################
############### Validacao do modelo PLS  ####################
#############################################################


validation_ggepls<-function(dados,config=NULL,genotype=NULL,env=NULL,resp=NULL, pred=NULL,c=2,tendency= T){
  
  if(is.null(config)== FALSE){
    genotype<- config$genotype
    env<- config$env
    resp<-config$resp
    pred<- config$pred
  }
 #dados<-df_cov
  na<-nrow(dados)- nrow(na.exclude(dados))
  if(na!=0){
    warning(paste0("\nDroping ",na," rows with NA") )
    dados<-na.exclude(dados)
  }
  names(dados)[match(genotype,names(dados))]<-"genotype"
  names(dados)[match(env,names(dados))]<-"env"
  gen<-as.character(unique(dados$genotype))
  dados$id<-1:nrow(dados)
  
  for (j in 1:length(gen)){
    
    if( j==1){
      r2.leave<-vector()
      ima.pred<-vector()
      predict<-vector()
      observ<-vector()
      geno<-vector()
      env<-vector()

    }
    
    cat("\n\nINICIANDO LEAVE ONE OUT:",gen[j],'\n')
    
    
    df.gen<-dados[dados$gen == gen[j],] #ACHAR UM JEITO DE AUTOMATIZAR
    predictors<-df.gen[pred]
    response <- df.gen[resp]
    
    total <- nrow(df.gen)
    # Barra de progresso
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    
    
    for( i in 1:total){
      
      pls <-plsreg1(predictors = predictors[-i,],response = response[-i,],comps = c,crosval = T)
      
      
      geno<-c(geno, gen[j])
      env<-c(env,df.gen[i,]$env)

      
      #Atualizando Barra de Progresso
      setTxtProgressBar(pb, i)
      
      pls.coefs<-pls$reg.coefs
      value<-t(c(1,as.numeric(predictors[i,])))%*%(pls$reg.coefs)
      predict<- c(predict,value)
      observ<- c(observ,response[i,])
      
    }
    
    if(j == length(gen)) df.loo<- data.frame(Genotype= geno, Env = env, Predicted = predict, Observed = observ, Residual = observ-predict)
  }
  
  require(dplyr)
  
  df.metrics<-df.loo %>% group_by(Genotype) %>% summarise(R2 = round(cor(Predicted, Observed) ^ 2,4),
                                                          RMSE = sqrt(mean((Predicted - Observed)^2)))
  
  
  if(tendency== T){
  (rmse.mean<-mean(df.metrics$RMSE))
  
  df.loo<-df.loo %>% mutate(Predictive_Tendency = case_when(Residual < -rmse.mean ~ "Overestimated",
                                                 Residual > rmse.mean ~ "Underestimated",
                                                 (Residual < rmse.mean) & ( Residual > - rmse.mean) ~ "Expected"))
  
  
  
  
  
  return(list(values = df.loo,metrics = df.metrics,rmse.mean= rmse.mean))
  
  }else{
    
    return(list(values = df.loo,metrics = df.metrics))
  }
}

#############################################################
############### Análise de covariáveis   ####################
#############################################################

covariates_analysis<-function(ggepls,rank,pred_names=NULL,resp_name=NULL){
  #rank=5;ggepls<-fit;resp_name = "MAI";pred_names=c("LATI","LONG","BDOD","CEC","CLAY","SAND",paste0("BIO",1:19),"ELEV")
  #rank=5;ggepls<-fit_sel;resp_name ="OI";pred_names= NULL
  
  if(ggepls$is.ggepls == T){
    
    gen<-names(ggepls$std_coef)
    
    #resp_name
    if(is.null(resp_name) ==T) resp_name <- ggepls$response
    
    #pred_names
    if(is.null(pred_names)== T) pred_names<- ggepls$predictors
    
    rownames(ggepls$std_coef)<-rownames(ggepls$coef)[-1]<-pred_names
    
    df.pred<-data.frame(var= c(ggepls$predictors,ggepls$response),
                        choosen_names = c(pred_names,resp_name))
    
    
    for(i in 1:length(gen)){
      
      if(i==1)g<-list()
      
      #Rankeando
      df.gen<-ggepls$std_coef[order(abs(ggepls$std_coef[,i]),decreasing = T),][gen[i]]
      best_cov<-rownames(df.gen)[1:rank]
      circle.gen<-ggepls$circle_coord[ggepls$circle_coord$gen == gen[i],]
      circle.gen<-merge(circle.gen,df.pred)
      
      #Tabela com as variáveis mais imoportantes
      df.rank<-data.frame(genotype = gen[i],rank=1:rank,covariates = best_cov,effects =ggepls$coef[best_cov,gen[i]], stdeffects = df.gen[best_cov,1])
      
      #cores
      circle.gen$highlight<-ifelse(circle.gen$choosen_names %in% c(best_cov,resp_name),circle.gen$choosen_names,"")
      circle.gen$covs<-ifelse(circle.gen$choosen_names %in% c(best_cov,resp_name),"",circle.gen$choosen_names)
      
      col.cov<-match(best_cov,circle.gen$highlight)
      col.resp<-match(resp_name,circle.gen$highlight)
      
      circle.gen$col<-"#0077b6";  circle.gen$col[col.cov]<-"black";circle.gen$col[col.resp]<-"orange"
      
      #Circulo de correlações
      g[[i]]<-ggplot(circle.gen,aes(x=t1,y=t2))+
        ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = 1),inherit.aes = FALSE,color="white")+
        geom_point(size=0.1)+
        geom_segment(aes(x = 0, y = 0, xend = t1, yend = t2), col = circle.gen$col)+
        labs(title = gen[i])+
        xlim(c(-1.2,1.2))+
        ylim(c(-1.2,1.2))+
        
        geom_text(aes(label=covs),size=6,col =circle.gen$col)+
        ggrepel::geom_label_repel(aes(label=highlight),size=6,col= circle.gen$col)+
        theme(plot.title = element_text(size = 26, face = "bold"),
              plot.subtitle =  element_text(size = 22),
              axis.title =element_text(size = 18),
              axis.text =element_text(size = 16),
              axis.text.x = element_blank(),
              panel.grid.minor = element_blank(),
              text = element_text(face="bold"))
      
      if(i==1) df.best<-df.rank else df.best<-rbind(df.best,df.rank)
    }
  }
  
  g[[1]]
  return(list(ranking = df.best,
              frequency = sort(table(df.best$covariates),decreasing = T),
              cicles = g))
}

######################################################################
########################  Corte de mapa  #############################
######################################################################

map_cut<-function(raster,shapefile){
  
  a<-crop(raster,shapefile)
  a<-mask(a,shapefile)
  
  return(a)
}

######################################################################
########################  Mapas temáticos  ###########################
######################################################################


ggepls_map<-function(ggepls, stack, cut = F,shapefile= NULL){
  
  #ggepls<-fit;stack<-e
  stack<-stack[[rownames(ggepls$coef)[-1]]] #Organizando na mesma ordem dos coeficientes
  gen<-ggepls$genotypes
  
  if(cut ==T){
    if(is.null(shapefile) ==F) stack<-map_cut(stack,shapefile) else warning("\nCorte nao realizado, insira o Shapefile")
  }
  
  for (j in 1:length(gen)){
    
    cat("\n\nPLOTAGEM DE MAPA:",as.character(gen[j]),'\n')
    
    pls.coefs<-ggepls$coef[,j] #Obtendo os coeficientes do genótipo
    
    
    total <- length(names(stack))
    # Barra de progresso
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    
    #plotagem do mapa
    for(i in 1:length(names(stack))){
      
      if(i==1) pred <-pls.coefs[1]+stack[[i]]*pls.coefs[i+1] else pred<-pred +stack[[i]]*pls.coefs[i+1]
      
      setTxtProgressBar(pb,i)
    }
    
    names(pred)<- gen[j]
    
    
    if(j==1) PRED<-pred else PRED<-stack(PRED,pred)
    
  }
  return(PRED)
}

######################################################################
########################  Quem vence onde  ###########################
######################################################################

get_X<-function(data,genotypes,region){
  #data<-data_cov; genotypes<-"MatGen_oficial"; env<-'unf'
  return(data.matrix(table(data[,genotypes],data[,region])>0))
}

rm_genotype<-function(raster,gen){
  for(i in 1:length(gen)) values(raster[[gen[i]]])<-NA
  
  return(raster)
}



which_won_where<- function(raster,rm_gen=NULL,by_region=F,shapefile=NULL,regions=NULL,X=NULL){
  
  #raster<-maps;rm_gen = c("CLZ003");regions="unf";by_region=T; X=X;shapefile=br
  #raster=maps
  #Removendo genótipos
  if(is.null(rm_gen)==F) raster<-rm_genotype(gen = rm_gen,raster = raster)
  
  
  #Fazendo quem vence onde para cara região
  
  df.gen<-data.frame(Genotype = names(raster), genotype_code =1:length(names(raster)) )
  
  if(by_region==T){
    
    reg<-shapefile@data[regions][,1]
    
    
    for( i in 1:length(reg)){
      sh.reg<-shapefile[shapefile@data[,regions] %in% reg[i],]
      gen.reg<-names(X[,reg[i]][X[,reg[i]]==T])
      
      if(length(setdiff(rownames(X),gen.reg))>0) raster.reg<-rm_genotype(raster,gen =  setdiff(rownames(X),gen.reg)) else raster.reg<-raster
      
      W<-whiches.max(raster.reg) #Quem ganha onde
      
      W<-map_cut(W,sh.reg)
      
      df.r<-data.frame(table(values(W))/sum(table(values(W))))
      names(df.r)[1]<-"genotype_code"
      df.r<-merge(df.r,df.gen)
      df.r$region<- reg[i]
      
      if(i==1){
        w<-W 
        df.reg<-df.r
        
      }else{
          
        w<-merge(w,W)
        df.reg<-rbind(df.reg,df.r)
        
        }
      
      #por regiao
      
      #total
      if(i==length(reg)){
        df.rec<-data.frame(table(values(w))/sum(table(values(w))))
        names(df.rec)[1]<-"genotype_code"
        
        df.rec<-merge(df.rec,df.gen)
        
        return(list(raster = w, general_recomendation = df.rec,region_recomendation = df.reg))
      }
      
    }
    
    
  }else{
    
    w<-which.max(raster)
    df.rec<-data.frame(table(values(w))/sum(table(values(w))))
    names(df.rec)[1]<-"genotype_code"
    
    df.rec<-merge(df.rec,df.gen)
    
    return(list(raster = w, general_recomendation = df.rec))
  }
}

#Plot com ggplot2
plot_recomendation<-function(w){
  
  w_pixel <- as(w$raster, "SpatialPixelsDataFrame")
  w_df <- as.data.frame(w_pixel)
  colnames(w_df) <- c("genotype_code", "x", "y")
  
  w_df<-merge(w_df,w$general_recomendation)
  
  g<-ggplot() +  
    geom_tile(data=w_df, aes(x=x, y=y, fill=Genotype))
  
  return(g)
}
