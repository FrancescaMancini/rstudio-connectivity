###script to extract lcm habitat from a polygon study area and calculate basic stats###
##Jenny Hodgson, 2/6/20

library(raster)
library(spdep)
library(rgeos)
library(rgdal)

setwd("/data/derived")
options(rasterTmpDir="../temprast")#not sure this has any effect
source("../notebooks/rstudio-connectivity/rstudio-connectivity/neighbwindowfuns.R")#load moving window functions for connectivity calculations

##########variables likely to change #########



#a shapefile with single polygon defining the study area
regionfname<- "derbyshire"  

#a raster with landcover codes, give name OR already have lcm loaded

#lcmfname<- "lcm07nenglandclip.tif"

lcm<- lcme

#working resolution

resn<- 25

#a sensible large square size to align to (round out to)

align<- 1000

#a short label for the focal habitat for naming files
hablabel<-"heath"
#what the focal habitat is called in the lookup table
habtext<- "heathland"

hlookup<- read.csv("../input-data/ToShare/Habitat_lookupTable.csv")

habcodes<- hlookup$LCM_class[hlookup$agg_class == habtext]


############preparing region data#############

#read in raster of all habitat
  #allhabitat<- raster(lcmfname)
  allhabitat<-lcm
#read in polygon of study area
bigpoly<- readOGR(dsn=getwd(), layer=regionfname)

gIsValid(bigpoly)
#[1] TRUE

#read in boundary as a line - may not be needed
#bigboundary<-gBoundary(bigpoly, byid=TRUE)


#set the sensible extent for all the rasters

algextent<- extent(bigpoly)

algextent@xmin<- floor(algextent@xmin/align)*align
algextent@ymin<- floor(algextent@ymin/align)*align
algextent@xmax<- ceiling(algextent@xmax/align)*align
algextent@ymax<- ceiling(algextent@ymax/align)*align


#create template empty raster

region1<- raster(algextent, resolution=resn, crs=crs(allhabitat))

#create raster mask of study area

regionp<- rasterize(bigpoly,region1) #slow

#regionp
#freq(regionp)
#     value   count
#[1,]     1 2084328
#[2,]    NA 1678872

regionb<- boundaries(regionp, type="outer") #slow

#freq(regionb)
#     value   count
#[1,]     0 2084328
#[2,]     1   10765
#[3,]   NaN 1668107

#alternative boundary - less useful as some cells may be both habitat and boundary
#regionbb<- rasterize(bigboundary,region1)

#freq(regionbb)
#     value   count
#[1,]     1   10764
#[2,]    NA 3752436



#freq(regionb*regionp)
#     value   count
#[1,]     0 2084328
#[2,]    NA 1678872


#crop habitat to the study region
habitatreg<- crop(allhabitat,algextent)

#freq(habitatreg)

#select the focal habitat

fhabitatreg<- habitatreg %in% habcodes
fhabitatreg<- mask(fhabitatreg, regionp) #mask after selecting habitat or NAs become 0s

#
freq(fhabitatreg)

#calculate vital statistics for habitat in region, units km

fhabitatbasics<-list(
  totarea=cellStats(fhabitatreg,sum)*resn*resn/1000/1000,
  proparea=cellStats(fhabitatreg,mean),
  regarea=cellStats(regionp,sum)*resn*resn/1000/1000,
  regperim=cellStats(regionb,sum)*resn/1000
)

fhabitatbasics$ccvc<- ccvc(raster=fhabitatreg,mask=regionp,resn=25,rescalc=200,pthresh=0.25,returnrasters=F)

mat<- expweightmat(dd=1000,resn=25,cutoff=0.975)

system.time({
  hk500<-hanski(raster=fhabitatreg,mask=regionp,dd=500,resn=25,cutoff=0.975,returnrasters=TRUE)
})

#user  system elapsed 
#426.184   0.000 426.155 

#plotting the Hanski connectivity surfaces - can comment these out to save time
plot(hk500$conr)
plot(hk500$conh)



system.time({
  hk2500<-hanski(raster=fhabitatreg,mask=regionp,dd=2500,resn=25,cutoff=0.975,returnrasters=TRUE)
})


#matrix 25* larger takes 16* longer

#add the Hanski overall metrics to the list of vital stats 
fhabitatbasics$hanski500<- hk500$coni
fhabitatbasics$hanski2500<- hk2500$coni


#save Rdata files that can be used as inputs to other analyses (e.g. adding habitat )
save (fhabitatbasics,fhabitatreg,regionp,regionb,
      regionfname,hablabel,habcodes,resn,align, file =paste(regionfname,hablabel,"Rdata",sep="."))
      
save (hk500,hk2500, file =paste(regionfname,hablabel,"hanski","Rdata",sep=".") )    
      