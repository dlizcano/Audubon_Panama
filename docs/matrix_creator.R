# function to make matrix per species per day, extracting and averaging hour 
# Modified from CI-TEAM Network code by Diego Lizcano
# October 2023
# Creates a nested list by species. One list presence another hour

f.matrix.creator2<-function(data,year){
  require(lubridate)
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
      mat_h[row,col]<-as.numeric((hours.cameras[j])) # remove hms
    }
    mat_h.nas<-is.na(mat)
    sum.nas<-apply(mat_h.nas,2,mean)# instead of sum
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

########################################


#code to shrink the matrix to exactly 9 columns: collapsing by 3 days
f.shrink.matrix.to9<-function(matrix){
  nc<-dim(matrix)[2]
  if(!nc%%9){ # of the number of columns is exactly divisible by 9
    newc<-nc%/%9
    old.cols<-seq(1,nc,newc)
    new.matrix<-matrix(NA,nr=nrow(matrix),nc=9)
    for(i in 1:9){
      new.matrix[,i]<-apply(matrix[,old.cols[i]:(old.cols[i]+newc-1)],1,max,na.rm=T)
    }
  } else{
    rem<-nc%%9
    newc<-nc%/%9
    old.cols<-seq(1,nc-rem,newc)
    new.matrix<-matrix(NA,nr=nrow(matrix),nc=9)
    for(i in 1:8)
      new.matrix[,i]<-apply(matrix[,old.cols[i]:(old.cols[i]+newc-1)],1,max,na.rm=T)
    new.matrix[,9]<-apply(matrix[,old.cols[9]:nc],1,max,na.rm=T) 
  }
  new.matrix[new.matrix=="-Inf"]<-NA
  rownames(new.matrix)<-rownames(matrix)
  new.matrix
}


#code to shrink the hour matrix to exactly 9 columns: collapsing by 3 days
f.shrink.matrix.h.to9<-function(matrix){
  nc<-dim(matrix)[2]
  if(!nc%%9){ # of the number of columns is exactly divisible by 9
    newc<-nc%/%9
    old.cols<-seq(1,nc,newc)
    new.matrix<-matrix(NA,nr=nrow(matrix),nc=9)
    for(i in 1:9){
      new.matrix[,i]<-apply(matrix[,old.cols[i]:(old.cols[i]+newc-1)],1,mean,na.rm=T)
    }
  } else{
    rem<-nc%%9
    newc<-nc%/%9
    old.cols<-seq(1,nc-rem,newc)
    new.matrix<-matrix(NA,nr=nrow(matrix),nc=9)
    for(i in 1:8)
      new.matrix[,i]<-apply(matrix[,old.cols[i]:(old.cols[i]+newc-1)],1,mean,na.rm=T)
    new.matrix[,9]<-apply(matrix[,old.cols[9]:nc],1,mean,na.rm=T) 
  }
  new.matrix[new.matrix=="-Inf"]<-NA
  rownames(new.matrix)<-rownames(matrix)
  new.matrix
}

