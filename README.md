# GIS-PLS - a pipeline to model the GxE interaction with Partial Least Square Regression and to predict and map the productivity in new areas with Geographic Information Systems

## Welcome to GISPLS!
The goal of this pipeline is to enhance the application of the GGEPLS approach to Geographic Information Systems (GIS) for predicting Genotype x Environment interaction (GxE). If you're interested, this repository contains the main script (GIS_PLS_en.R or GIS_PLS_ptbr.R) with all the necessary instructions to conduct the analysis using inventory data (production sites). The pipeline functions were designed to make the approach user-friendly, even for those without a strong background in GIS.

## Inventory data
The provided data shows standardized predictions of the Mean Annual Increment (MAI) of wood at 7 years old, based on productive eucalyptus sites across major regions of Brazil.

## Environmental data
The GISPLS approach is a versatile tool that can handle various types of environmental data. In our pipeline, there are two functions, "get_WC" and "get_SG", which can quickly provide rasters from WorldClim and SoilGrids at four different spatial resolutions: 10 minutes, 5 minutes, 2.5 minutes, and 30 seconds. Please extract the files from "Functions_GISPLS.zip" and move them to the same location as the main script.


## Functions
### covariates_analysis
covariates_analysis(ggepls,rank,pred_names=NULL,resp_name=NULL)

##### ggepls: is obtained by the ggepls_fit() function.
##### rank: number of covariates that will be ranked.
##### pred_names: covariates personalized names.
##### resp_names: response variable name.

### fit_ggepls
fit_ggepls (data,config=NULL,genotype=NULL, env=NULL,resp=NULL, pred=NULL,c=2)

##### data: data frame containing the genotype and environmental information, such as response and covariates indexes.
##### config: shortcut list to set genotype, env, resp, and pred.
##### genotype: genotype column name.
##### env: environment or regions column name.
##### resp: response variable index.
##### pred: covariates indexes.

### get_X
get_X(data,genotypes,region)

##### data: data frame containing genotypes and region information.
##### genotype: genotype column name.
##### region: environment or regions column name.

### getLonLat 
getLonLat(raster)
  
##### raster: template raster to extract the average longitude and latitude
### getSG
getSG (attributes = NULL, layers ="30-60cm_mean", resample = F, raster=NULL,tr=c(4500,4500)

##### attributes: soil attributes that should be downloaded. Options: bdod, cec, cfvo, clay, nitrogen, ocd, ocs, phh2o, sand, silt, soc, wrb. For more information access: https://www.isric.org/explore/soilgrids/faq-soilgrids
##### layers: soil layers that sould be downloaded. Options: 0-5cm_mean, 5-15cm_mean, 15-30cm_mean, 30-60cm_mean, 60-100cm_mean, 100-200cm_mean. It's also possible to download other quantiles, for example, 0-5cm_Q0.05, 0-5cm_Q0.5, 0-5cm_Q0.95, 0-5cm_uncertainty. For more information access: https://www.isric.org/explore/soilgrids/faq-soilgrids
##### resample: Soilgrid rasters are composed by the Homolosine projection. If TRUE you will need to provide a template raster with the desired DATUM, extended, and special resolution.
##### raster: template raster.
##### tr: special resolution in meters required by rgdal to download Soilgrids rasters

### getWC
getWC(var=NULL,res=NULL)

##### var: climate variables that should be downloaded. Options: tmin, tmax, tavg, prec, srad, wind, vapr, bioc, elev.
##### res: the spacial resolution that rasters should be downloaded. Options: r10min, r5min, r2.5min, r30sec


### ggepls_map
ggepls_map(ggepls, stack, cut = F,shapefile= NULL)

##### ggepls: ggepls object obtained by ggepls_fit() function.
##### stack: stack raster containing all covariates rasters
##### cut: If TRUE, you will need to provide a shapefile template to cut the rasters.
##### shapefile: template shapefile

### grid_point
grid_point (raster,point)
##### raster: template raster where the grid will be obtained
##### point: data points from your data

### map_cut 
  map_cut(raster,shapefile)

This function is simply a combination of the crop and mask functions from the raster package.
##### raster: raster to be cut
##### shapefile: template shapefile (or raster)

### plot_recomendation
plot_recomendation(w)

##### w: object obtained from which_won_where() function

### rastertype
rastertype( stack,points,point_id=NULL)

##### stack: stack raster containing all covariates rasters
##### points: data points from your data
##### point_id: if NULL an extra column with an id to each point will be created


### rm_genotype
rm_genotype(raster,gen)

##### raster: raster containing all mapped genotype's performance
##### gen: genotype’s name to be removed

### stack_rasters
stack_rasters (path,pattern)

This function is simply a combination of the list.files() function from the base package and the raster() function from the raster package.

##### path: rasters files path
##### pattern: file name pattern

### validation_ggepls
validation_ggepls (data,config=NULL,genotype=NULL,env=NULL,resp=NULL, pred=NULL,c=2,tendency= T)

##### data: data frame containing the genotype and environmental information, such as response and covariates indexes
##### config: shortcut list to set genotype, env, resp, and pred
##### genotype: genotype column name
##### env: environment or regions column name
##### resp: response variable index
##### pred: covariates indexes
##### c: number of latent variables
##### tendency: if TRUE an extra column with an error classification will be created

### which_won_where
  which_won_where(raster,rm_gen=NULL,by_region=F,shapefile=NULL,regions=NULL,X=NULL)

##### raster: raster containing all mapped genotype's performance
##### rm_gen: genotype’s name to be removed
##### by_region: if TRUE the X matrix from get_X() function will be used to compare only genotypes tested by region
##### shapefile: shapefile that contains region information
##### regions: attribute name in shapefile that identifies the regions
##### X: matrix from get_X()










