##########################################################
######## Grid para obter as medias por quadrante #########
##########################################################

grid_point<-function(raster,point,point_id){
  #point<-WORLD
  #raster<- raster("./Rasters_2.5m/BR_wc2.1_2.5m_elev.tif")
  #point_id<-"id"
  
  values(raster)<-1:length(raster)
  point@data$id<-point@data[,point_id]
  
  #Barra de progresso
  total <- length(point@data$id)
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
  #stack = e; points =inv[1:5,]
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
  
  return(edacli)
}

##########################################################
############### Ajuste do modelo PLS  ####################
##########################################################

fit_ggepls<-function(dados,config=NULL,genotype=NULL,env=NULL,resp=NULL, pred=NULL,c=2){
  
  if(is.null(config)== FALSE){
    genotype<- config$genotype
    env<- config$env
    resp<-config$resp
    pred<- config$pred
  }
  
  #dados<-data_cov;genotype <- "MatGen_oficial";env <- "unf";pred <- c(6:7,14:38); resp<- 5 
  names(dados)[match(genotype,names(dados))]<-"genotype"
  names(dados)[match(env,names(dados))]<-"env"
  gen<-unique(dados$genotype)
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
    
    cat("\nINICIANDO AJUSTE GGEPLS:",gen[j],'\n')
    
    df.gen<-dados[dados$gen == gen[j],] 
    predictors<-df.gen[pred]
    response <- df.gen[resp]
    
    total <- nrow(df.gen)
    # Barra de progresso
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    
    #PLS
    pls <-plsdepot::plsreg1(predictors = predictors,response = response,comps = c,crosval = T)
    
    #Latent variables
    coo<-as.data.frame(pls$cor.xyt[,c('t1','t2')])
    coo$var<-rownames(coo)
    coo$gen <-gen[j]
    circle <- coo
    
    geno<-c(geno, df.gen$genotype)
    env<-c(env,df.gen$env)
    point<-c(point,df.gen$id)
    
    #Atualizando Barra de Progresso
    setTxtProgressBar(pb, total)
    
    pls.coefs<-pls$reg.coefs
    value<-as.numeric(data.matrix(cbind(1,predictors))%*%(pls$reg.coefs))
    predict<- c(predict,value)
    observ<- c(observ,response[,1])
    
    if(j==1) {
      beta<-pls.coefs
      beta_std<-pls$std.coefs
      CIRCLE<-circle
      }else
        {
      beta<-cbind(beta,pls.coefs)
      beta_std<- cbind(beta_std,pls$std.coefs)
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


validation_ggepls<-function(dados,config=NULL,genotype=NULL,env=NULL,resp=NULL, pred=NULL,c=2){
  
  #dados<-data_cov;genotype <- "MatGen_oficial";env <- "unf";pred <- c(6:7,14:38); resp<- 5 
  if(is.null(config)== FALSE){
    genotype<- config$genotype
    env<- config$env
    resp<-config$resp
    pred<- config$pred
  }
  
  names(dados)[match(genotype,names(dados))]<-"genotype"
  names(dados)[match(env,names(dados))]<-"env"
  gen<-unique(dados$genotype)
  dados$id<-1:nrow(dados)
  
  for (j in 1:length(gen)){
    
    if( j==1){
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
    
    cat("\nINICIANDO LEAVE ONE OUT:",gen[j],'\n')
    
    
    df.gen<-dados[dados$gen == gen[j],] #ACHAR UM JEITO DE AUTOMATIZAR
    predictors<-df.gen[pred]
    response <- df.gen[resp]
    
    total <- nrow(df.gen)
    # Barra de progresso
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    
    
    for( i in 1:total){
      
      pls <-plsdepot::plsreg1(predictors = predictors[-i,],response = response[-i,],comps = c,crosval = T)
      
      
      geno<-c(geno, gen[j])
      env<-c(env,df.gen[i,]$env)
      point<-c(point,df.gen[i,]$id)
      
      #Atualizando Barra de Progresso
      setTxtProgressBar(pb, i)
      
      pls.coefs<-pls$reg.coefs
      value<-t(c(1,as.numeric(predictors[i,])))%*%(pls$reg.coefs)
      predict<- c(predict,value)
      observ<- c(observ,response[i,])
      
      
    }
    
    
    if(j == length(gen)) df.loo<- data.frame(Genotype= geno, Env = env, Predicted = predict, Observed = observ)
  }
  
  require(dplyr)
  
  df.metrics<-df.loo %>% group_by(Genotype) %>% summarise(R2 = round(cor(Predicted, Observed) ^ 2,4),
                                                          RMSE = sqrt(mean((Predicted - Observed)^2)))
  return(list(values = df.loo,metrics = df.metrics))
}

#############################################################
############### Análise de covariáveis   ####################
#############################################################
library(plsdepot)
library(plotrix)
library(ggplot2)


covariates_analysis<-function(ggepls,rank,pred_names=NULL,resp_name=NULL){
  #rank=5;ggepls<-fit;resp_name = "MAI";pred_names=c("LATI","LONG","BDOD","CEC","CLAY","SAND",paste0("BIO",1:19),"ELEV")
  #rank=5;ggepls<-fit;resp_name ="OI";pred_names= NULL
  if(ggepls$is.ggepls == T){
    
    gen<-names(ggepls$std_coef)
    
    #resp_name
    if(is.null(resp_name) ==T) resp_name <- ggepls$response
    
    
    #pred_names
    if(is.null(pred_names)== T) pred_names<- ggepls$predictors
    
    rownames(ggepls$std_coef)<-rownames(ggepls$coef)[-1]<-pred_names
    ggepls$circle_coord$var<-c(pred_names,resp_name)

    for(i in 1:length(gen)){
      if(i==1)g<-list()

      
      #Rankeando
      df.gen<-ggepls$std_coef[order(abs(ggepls$std_coef[,i]),decreasing = T),][gen[i]]
      best_cov<-rownames(df.gen)[1:rank]
      circle.gen<-ggepls$circle_coord[ggepls$circle_coord$gen == gen[i],]
      
      #Tabela com as variáveis mais imoportantes
      df.rank<-data.frame(genotype = gen[i],rank=1:rank,covariates = best_cov,effects =ggepls$coef[best_cov,gen[i]], stdeffects = df.gen[best_cov,1])
      

      #Renomeando as variáveis no circulo
    
      
      circle.gen$rank<-ifelse(circle.gen$var %in% c(best_cov,resp_name),circle.gen$var,"")
      circle.gen$unrank<-ifelse(circle.gen$var %in% c(best_cov,resp_name),"",circle.gen$var)
      col.end<-match(best_cov,circle.gen$rank)
      
      circle.gen$col<-"#0077b6";  circle.gen$col[col.end]<-"black";circle.gen$col[nrow(circle.gen)]<-"orange"
      
      #Circulo de correlações
      g[[i]]<- ggplot(circle.gen,aes(x=t1,y=t2))+
        geom_circle(aes(x0 = 0, y0 = 0, r = 1),inherit.aes = FALSE,color="white")+
        geom_point(size=0.1)+
        geom_segment(aes(x = 0, y = 0, xend = t1, yend = t2), col = circle.gen$col)+
        labs(title = gen[i])+
        xlim(c(-1.2,1.2))+
        ylim(c(-1.2,1.2))+
        
        geom_text(aes(label=unrank),size=6,col =circle.gen$col)+
        geom_label_repel(aes(label=rank),size=6,col= circle.gen$col)+
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

  return(list(ranking = df.best, cicles = g))
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


ggepls_map<-function(ggepls, stack, cut = F,shapefile= NULL ){
  
  #ggepls<-fit;stack<-e
  stack<-stack[[rownames(ggepls$coef)[-1]]] #Organizando na mesma ordem dos coeficientes
  gen<-ggepls$genotypes

  
  for (j in 1:length(gen)){
    
    cat("\n\nPLOTAGEM DE MAPA:",gen[j],'\n')
    
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
    

    if(cut ==T){
      if(is.null(shapefile) ==F) pred<-map_cut(pred,shapefile) else warning("\nCorte nao realizado, insira o Shapefile")
    }
    
    if(j==1) PRED<-pred else PRED<-stack(PRED,pred)
    
  }
  return(PRED)
}

######################################################################
########################  Quem vence onde  ###########################
######################################################################

get_X<-function(data,genotypes,region){
  #data<-data_cov; genotypes<-"MatGen_oficial"; env<-'unf'
  return(data.matrix(table(data_cov[,genotypes],data_cov[,region])>0))
}

rm_genotype<-function(raster,gen){
    for(i in 1:length(gen)) values(raster[[gen[i]]])<-NA
    
    return(raster)
}
  


which_won_where<- function(raster,rm_gen=NULL,by_region=F,shapefile=NULL,regions=NULL,X=NULL){
  
  #raster<-maps;rm_gen = c("CLZ003");regions="unf";by_region=T
  
  #Removendo genótipos
  if(is.null(rm_gen)==F) raster<-rm_genotype(gen = rm_gen,raster = raster)
  
  
  #Fazendo quem vence onde para cara região
  
  
  if(by_region==T){
    
    reg<-br@data[regions][,1]
    
    for( i in 1:length(reg)){
      
      sh.reg<-br[br@data[,regions] %in% reg[i],]
      gen.reg<-names(X[,reg[i]][X[,reg[i]]==T])
      
      if(length(setdiff(rownames(X),gen.reg))>0) raster.reg<-rm_genotype(raster,gen =  setdiff(rownames(X),gen.reg)) else raster.reg<-raster
      
      W<-whiches.max(raster.reg) #Quem ganha onde
      
      W<-map_cut(W,sh.reg)
      
      if(i==1)w<-W else w<-merge(w,W)
    }
    
    
  }else{
    
    w<-which.max(raster)
  }
  return(w)
  
}




