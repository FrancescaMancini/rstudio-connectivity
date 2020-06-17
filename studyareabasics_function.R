## Function for study area basic stats, hanski metric 
## and climate change vulnerability tool connectivity metric

studyareabasics <- function(inputs){
  
  # Tell cluster where your libraries are
  
  .libPaths(c("/data/conda/myRpy/lib/R/library",.libPaths()))
  Sys.setenv(GDAL_DATA = "/data/conda/myRpy/share/gdal")
  Sys.setenv(PROJ_LIB = "/data/conda/myRpy/share/proj")
  
  
  # Set cluster working directory to source notebook
  setwd(inputs$mywd)
  
  # load libraries
  library(raster)
  library(spdep)
  library(rgeos)
  library(rgdal)
  
  #load moving window functions for connectivity calculations
  source("neighbwindowfuns.R")
  
  sink(file.path(getwd(), paste0('console_',inputs$regionfname,'_', inputs$hablabel,'.txt')))
  
  
  hlookup<- read.csv(file = inputs$hlookuppath)
  
  if(inputs$habtext == "semi-natural"){
    habcodes <- hlookup$LCM_class[hlookup$semi_natural == "YES"]
  }
  else{habcodes<- hlookup$LCM_class[hlookup$agg_class == inputs$habtext]}
  
  
  
  ############preparing region data#############
  
  # read in LCM
  
  lcm <- raster(inputs$rasterpath)
  
  cat(paste0("\n", 'LCM loaded', "\n", Sys.time(), "\n"))
  
  #read in polygon of study area
  bigpoly<- spTransform(readOGR(dsn= inputs$regionspath, layer=inputs$regionfname), crs(lcm))
  
  cat(paste0("\n", 'Bigpoly loaded', "\n", Sys.time(), "\n"))
  
  #set the sensible extent for all the rasters
  
  algextent<- extent(bigpoly)
  
  algextent@xmin<- floor(algextent@xmin/inputs$align)*inputs$align
  algextent@ymin<- floor(algextent@ymin/inputs$align)*inputs$align
  algextent@xmax<- ceiling(algextent@xmax/inputs$align)*inputs$align
  algextent@ymax<- ceiling(algextent@ymax/inputs$align)*inputs$align
  
  cat(paste0("\n", 'Extent', "\n", Sys.time(), "\n"))
  
  #create template empty raster
  
  region1<- raster(algextent, resolution=inputs$resn, crs=crs(lcm))
  
  cat(paste0("\n", 'Raster created', "\n", Sys.time(), "\n"))
  #create raster mask of study area
  
  regionp<- rasterize(bigpoly,region1) #slow
  regionb<- boundaries(regionp, type="outer") #slow
  
  cat(paste0("\n", 'Mask created', "\n", Sys.time(), "\n"))
  
  #crop habitat to the study region
  habitatreg<- crop(lcm,algextent)
  
  cat(paste0("\n", 'Habitat cropped', "\n", Sys.time(), "\n"))
  
  #select the focal habitat
  
  fhabitatreg<- habitatreg %in% habcodes
  fhabitatreg<- mask(fhabitatreg, regionp) #mask after selecting habitat or NAs become 0s
  
  cat(paste0("\n", 'Focal habitat', "\n", Sys.time(), "\n"))
  
  
  #calculate vital statistics for habitat in region, units km
  
  fhabitatbasics<-list(
    totarea=cellStats(fhabitatreg,sum)*inputs$resn*inputs$resn/1000/1000,
    proparea=cellStats(fhabitatreg,mean),
    regarea=cellStats(regionp,sum)*inputs$resn*inputs$resn/1000/1000,
    regperim=cellStats(regionb,sum)*inputs$resn/1000
  )
  
  cat(paste0("\n", 'Stats calculated', "\n", Sys.time(), "\n"))
  
  # connectivity metric from Climate Change vulnerability tool
  
  fhabitatbasics$ccvc<- ccvc(raster=fhabitatreg,mask=regionp,resn=25,rescalc=200,pthresh=0.25,returnrasters=F)
  
  cat(paste0("\n", 'ccvc calculated', "\n", Sys.time(), "\n"))
  
  # hanski metric
  
  hk500<-hanski(raster=fhabitatreg,mask=regionp,dd=500,resn=25,cutoff=0.975,returnrasters=TRUE)
  
  cat(paste0("\n", 'Hanski500 computed', "\n", Sys.time(), "\n"))
  
  
  hk2500<-hanski(raster=fhabitatreg,mask=regionp,dd=2500,resn=25,cutoff=0.975,returnrasters=TRUE)
  
  cat(paste0("\n", 'Hanski2500 computed', "\n", Sys.time(), "\n"))
  
  #add the Hanski overall metrics to the list of vital stats 
  fhabitatbasics$hanski500<- hk500$coni
  fhabitatbasics$hanski2500<- hk2500$coni
  
  # save files
  
  save (fhabitatbasics,fhabitatreg,regionp,regionb,
        inputs$regionfname,inputs$hablabel,habcodes,
        resn,align, file =paste0("/data/derived/", inputs$regionfname, "/",
                                 paste(inputs$regionfname,inputs$hablabel,"Rdata",sep=".")))
  
  save (hk500,hk2500, file =paste0("/data/derived/", inputs$regionfname, "/", 
                                   paste(inputs$regionfname,inputs$hablabe,"hanski","Rdata",sep=".")))    
  
  # return "job done" msg and date and time and close the connection
  cat(paste0("\n", 'Job done', "\n", Sys.time(), "\n"))
  
  sink()
  
}