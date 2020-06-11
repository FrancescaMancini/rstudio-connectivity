#####functions for moving window landscape evaluation#####
##Jenny Hodgson, May 2020

hanski<- function(raster,mask,dd,resn,cutoff=0.975,returnrasters=TRUE){
  
  mat<- expweightmat(dd=dd,resn=resn,cutoff=cutoff)
  
  conr<- focal(x=raster, w=mat,  na.rm=TRUE)
  
  #only count connectivity of habitat cells
  conh<- conr*raster
  
  #average connectivity of habitat, where 1 would be found if 100% coverage
  coni<- cellStats(conh,mean)
  
  if(returnrasters){
    return(list(conr=conr,conh=conh,coni=coni))
  }else{
    return(coni)
  }
}

expweightmat<- function(dd,resn,cutoff=0.975){
  #all distances should be in metres, if to be used with a raster in m
  alpha<- 2/dd
  
  #find correct cutoff distance, as alpha*r
  
  rem<- 1-cutoff
  
  try<-seq(-log(rem),pmax(3,-log(rem)*2), -log(0.95))
  
  propk<- exp(-try)*(try+1)
  
  rdd<- try[order(abs(propk-rem))[1]]
  
  #radius in m
  rd<- rdd/alpha
  
  #radius in no. of cells
  rc<- round(rd/resn)
  
  mat<- sqrt(outer((0:rc)^2, (0:rc)^2, FUN = "+"))
  
  fact<- resn*resn* alpha^2/2/pi
  mat<- fact*exp(-mat*resn*alpha)
  
  mat[mat< (fact*exp(-rdd)) ]<- 0
  #don't count self
  mat[1,1]<-0
  
  #construct 4 quadrants
  mat2<- cbind( mat[,(rc+1):2],mat)
  mat4<- rbind( mat2[(rc+1):2,],mat2)
  
  #rescale to sum to 1
  mat4<- mat4/sum(mat4)
  
  mat4}


ccvc<- function(raster,mask,resn,rescalc=200,pthresh=0.25,returnrasters=TRUE){
  # Note mask is not used yet - we rely on fact that raster is NA outside mask
  
  # Create matrix show in Taylor,Knight and Harfoot 2014, Figure5
  matq<- matrix(c(0,6,3,6,4,2,3,2,1),nrow=3)
  mat2<- cbind( matq[,3:2],matq)
  mat<- rbind( mat2[3:2,],mat2)
  #rescale to sum to 1
  mat<- mat/sum(mat)
  
  #aggregation of raster
  aggfact<-round(rescalc/resn)
  coarseraster<- aggregate(raster,fact=aggfact,fun=sum,na.rm=TRUE)
  #convert back to pres/absence using a threshold 
  coarseraster<- (coarseraster/aggfact/aggfact) > pthresh
  
  #do moving window
  conr<- focal(x=coarseraster, w=mat,  na.rm=TRUE)
  
  #only count connectivity of habitat cells
  conh<- conr*coarseraster
  
  #average connectivity of habitat, where 1 would be found if 100% coverage
  coni<- cellStats(conh,mean)
  
  if(returnrasters){
    return(list(conr=conr,conh=conh,coni=coni))
  }else{
    return(coni)
  }
}	