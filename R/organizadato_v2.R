

f.matrix.creator2<-function(data,year){
  #results object
  res<-list()
  
  #get the dimensions of the matrix
  
  #list if sanpling units
  cams<-unique(data$locationID)
  cams<-sort(cams)
  rows<-length(cams)
  species<-unique(data$scientificName)
  #start and end dates of sampling periods
  # data<-data[data$Sampling.Period==year,]
  min<-min(as.Date(as.character(data$start_date), "%Y-%m-%d"))
  max<-max(as.Date(as.character(data$end_date), "%Y-%m-%d"))
  cols<-max-min+1
  
  #sampling period
  date.header<-seq(from=min,to=max, by="days")
  
  mat<-matrix(NA,rows,cols,dimnames=list(cams,as.character(date.header)))
 
  #for all cameras, determine the open and close date and mark in the matrix
  start.dates<-tapply(as.character(data$start_date),data$locationID,unique)
  nms<-names(start.dates)
  # start.dates<-ymd(start.dates)
  names(start.dates)<-nms
  end.dates<-tapply(as.character(data$end_date),data$locationID,unique)
  # end.dates<-ymd(end.dates)
  names(end.dates)<-nms
  
  #outline the sampling periods for each camera j
  for(j in 1:length(start.dates)){
    #for each camera beginning and end of sampling
    low<-which(date.header==as.Date(as.character(start.dates[j]), format = "%Y-%m-%d"))
    hi<-which(date.header==as.Date(as.character(end.dates[j]), format = "%Y-%m-%d"))
    if(length(low)+length(hi)>0){
      indx<-seq(from=low,to=hi)
      mat[names(start.dates)[j],indx]<- 0
    } else next
  }
  mat.template<-mat
  mat_h <- mat
  #get the species
  #species<-unique(data$bin)
  #construct the matrix for each species i
  for(i in 1:length(species)){
    indx<-which(data$scientificName==species[i])
    #dates and cameras when/where the species was photographed
    dates<-data$eventDate[indx]
    cameras<-data$locationID[indx]
    hora <- data$time[indx]
    
    dates.cameras <- data.frame(dates, cameras)
    hours.cameras <- data.frame(hora, cameras)
    #unique combination of dates and cameras 
    dates.cameras<-unique(dates.cameras)
    hours.cameras <- unique(hora, cameras)
    #fill in the obs matrix
    for(j in 1:length(dates.cameras[,1])){
      col<-which(date.header==as.character( dates.cameras[j,1]))
      row<-which(cams==as.character( dates.cameras[j,2]))
      mat[row,col]<-1
    }
    mat.nas<-is.na(mat)
    sum.nas<-apply(mat.nas,2,sum)
    indx.nas<-which(sum.nas==rows)
    if(length(indx.nas)>0){
      mat<-mat[,-indx.nas]
    }
    
    
    #fill in the hour matrix
    for(j in 1:length(dates.cameras[,1])){
      col<-which(date.header==as.character( dates.cameras[j,1]))
      row<-which(cams==as.character( dates.cameras[j,2]))
      mat_h[row,col]<-as.numeric(hms(hours.cameras[j]))
    }
    mat_h.nas<-is.na(mat)
    sum.nas<-apply(mat_h.nas,2,sum)
    indx.nas<-which(sum.nas==rows)
    if(length(indx.nas)>0){
      mat_h<-mat_h[,-indx.nas]
    }

    occur<-c(res,list(mat)) #lista anidada
    dete <- c(res,list(mat_h))
    res<-c(res,list(mat))
    # res<-list(occur=list(mat), hora=list(mat_h), sp=list(species)) #lista anidada
    #return the matrix to its original form
    mat<-mat.template
    
  }
  
  # names(occur)<-species
  names(dete)<-species
  #res<-lapply(res,f.dum)
  return(list(occur=occur, dete=dete))
  
}




#code to shrink the matrix to exactly 26 columns: Aprox una semana
f.shrink.matrix.to26<-function(matrix){
  nc<-dim(matrix)[2]
  if(!nc%%26){ # of the number of columns is exactly divisible by 15
    newc<-nc%/%26
    old.cols<-seq(1,nc,newc)
    new.matrix<-matrix(NA,nr=nrow(matrix),nc=26)
    for(i in 1:26){
      new.matrix[,i]<-apply(matrix[,old.cols[i]:(old.cols[i]+newc-1)],1,max,na.rm=T)
    }
  } else{
    rem<-nc%%26
    newc<-nc%/%26
    old.cols<-seq(1,nc-rem,newc)
    new.matrix<-matrix(NA,nr=nrow(matrix),nc=26)
    for(i in 1:25)
      new.matrix[,i]<-apply(matrix[,old.cols[i]:(old.cols[i]+newc-1)],1,max,na.rm=T)
    new.matrix[,26]<-apply(matrix[,old.cols[26]:nc],1,max,na.rm=T) 
  }
  new.matrix[new.matrix=="-Inf"]<-NA
  rownames(new.matrix)<-rownames(matrix)
  new.matrix
}




############### load data

library(readr)
library(tidyverse)
parita_data_full <- read_delim("C:/CodigoR/AudubonPanama/data/all_records_audubon_panama_to_occu.csv", 
                       delim = "\t", escape_double = FALSE, 
                       col_types = cols(eventDate = col_date(format = "%Y-%m-%d")), 
                       trim_ws = TRUE)


# filter by threshold
parita_data <-  parita_data_full %>% filter(confidence >= threshold)


# insert Max and min 
# start_date
parita_data2 <- parita_data %>% 
  group_by(locationID) %>% 
  mutate(start_date = min(eventDate),
         end_date= max(eventDate))
                          
parita <- f.matrix.creator2(data=parita_data2, year=2023)
covs <- read_csv("C:/CodigoR/AudubonPanama/shp/sites_covs_parita_nona.csv")
library(unmarked)
library(sf)
library(mapview)

umf_y_full<- unmarkedFrameOccu(y= parita$dete[[5]]) #Actitis macularius
siteCovs(umf_y_full) <- data.frame(elevation =  scale(as.numeric(covs$Altitude)),
                                   AGB_Spawn = scale(as.numeric(covs$AGB_Spawn)),
                                   BGB_Spawn = scale(as.numeric(covs$BGB_Spawn)),
                                   roads = scale(as.numeric(covs$roads_final_v2)),
                                   carbon = scale(as.numeric(covs$organic_carbon_stock_0.30m)),
                                   NDVI=scale(covs$S2_NDVI_median_v2)
                                   ) # data.frame(Elev=full_covs$Elev) # Full

covariables <- data.frame(elevation =  scale(as.numeric(covs$Altitude)),
                          AGB_Spawn = scale(as.numeric(covs$AGB_Spawn)),
                          BGB_Spawn = scale(as.numeric(covs$BGB_Spawn)),
                          roads = scale(as.numeric(covs$roads_final_v2)),
                          carbon = scale(as.numeric(covs$organic_carbon_stock_0.30m)),
                          NDVI=scale(covs$S2_NDVI_median_v2),
                          Station=covs$Name,
                          Lon=covs$Longitude,
                          Lat=covs$Latitude)


#effort <- list(hour= as.data.frame(Camptostoma_obsoletum_occu_hora[,10:36]))

# obsCovs(umf_y_full) <- list(hour= as.data.frame(Camptostoma_obsoletum_occu_hora[,10:36]))

#######Graficar umf
plot(umf_y_full) 



#########################
library(camtrapR)
library(purrr)
library(DT)

# data("camtraps")


#Replace Nas by 0 - Another way
for(j in 1:length(parita$dete)){
  parita$dete[[j]] <- scale(parita$dete[[j]] %>% replace(is.na(.), 0))
}


effort <- matrix(1, nrow = nrow(parita$occur[[1]]), ncol =  ncol(parita$occur[[1]])) # no effect

# bundle data
data_list <- list(ylist    = parita$occur,
                  siteCovs = covariables,
                  obsCovs  = list(hour = parita$dete$`Zenaida asiatica`,
                                  effort = effort)
                  )
                  # )  # is identical for all species



# text file to save the model
modelfile1 <- tempfile(fileext  = ".txt")


model1_jags <- communityModel(
  occuCovs = list(ranef = c("NDVI", "AGB_Spawn", "roads")), 
  detCovsObservation = list(ranef = "hour"),
  intercepts = list(det = "fixed", occu = "ranef"),
  data_list = data_list,
  modelFile = modelfile1)

# model2_jags <- communityModel(
#   occuCovs = list(ranef = c("NDVI", "AGB_Spawn", "roads")), 
#   detCovs = list(ranef = "hour"),
#   intercepts = list(det = "fixed", occu = "ranef"),
#   data_list = data_list,
#   modelFile = modelfile2,
#   effortCov = "effort")

library(rjags)
library(nimble)

### run
out_ahm_jags <- fit(model1_jags, 
                    n.iter = 100000, 
                    n.burnin = 5000,
                    thin = 5,
                    chains = 3,
                    quiet = T
)

# plot
plot_resp_jags_occu <- plot_effects(model1_jags, 
                                    out_ahm_jags)
plot_resp_jags_occu

# coefficients
plot_eff_jags_occu <- plot_coef(model1_jags, 
                                out_ahm_jags)
plot_eff_jags_occu


##### detection
plot_coef(model1_jags, 
          out_ahm_jags, 
          submodel = "det")   

#####################################################
# Save an object to a file
saveRDS(out_ahm_jags, file = "out_ahm_jags.rds")
# Restore the object
readRDS(file = "out_ahm_jags.rds")

####################################################
#### hora cuadratica
data_list$obsCovs$hour_squared <- data_list$obsCovs$hour ^2

model2_jags <- communityModel(
  occuCovs = list(ranef = c("NDVI", "AGB_Spawn", "roads")), 
  detCovsObservation = list(ranef = c("hour", "hour_squared")),
  intercepts = list(det = "ranef", occu = "ranef"),
  data_list = data_list,
  effortCov = "effort",
  nimble = FALSE,
  #keyword_quadratic = "_squared",
  modelFile = modelfile1)


### run
out_ahm_jags2  <- fit(model2_jags, 
                    n.iter = 100000, 
                    n.burnin = 10000,
                    thin = 5,
                    chains = 3,
                    quiet = T,
                    cores = 4
)

# traceplots
plot(fit.nimble.comp)
 
# plot response curves (= marginal effect plots)
plot_resp_jags_occu <- plot_effects(model2_jags, 
                                      out_ahm_jags2)
plot_resp_jags_occu

plot_resp_jags_det <-plot_effects(model2_jags, 
             out_ahm_jags2, submodel = "det")

plot_resp_jags_det

# coefficients
plot_eff_jags_occu <- plot_coef(model2_jags, 
                                  out_ahm_jags2)
plot_eff_jags_occu


##### detection
plot_coef(model2_jags, 
          out_ahm_jags2, 
          submodel = "det")   

