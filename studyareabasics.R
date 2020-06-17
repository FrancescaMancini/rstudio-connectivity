###script to extract lcm habitat from a polygon study area and calculate basic stats###
##Jenny Hodgson, 2/6/20
## Updated 17/06/2020

## this was run on a spark cluster 

rm(list=ls())

library(SparkR, lib.loc = file.path(Sys.getenv('SPARK_HOME'), "R", "lib"))

# prepare input parameters
# needs to be a list of lists
sparkInputs <- lapply(1:12, 
                      FUN = function(i){list(index = i, regionfname = rep(c("derbyshire", "w_yorkshire", "w_sussex"), each = 4)[i],
                                             regionspath = "/data/derived",
                                             rasterpath = "/data/input-data/ToShare/LCM2015_England_cropped/LCM2015_eng.tif",
                                             resn = 25, align = 1000, 
                                             hablabel = rep(c("heath", "wood", "grass", "snat"), 3)[i],
                                             habtext = rep(c("heathland", "woodland", "semi-natural grassland", "semi-natural"), 3)[i], 
                                             hlookuppath = "/data/input-data/ToShare/Habitat_lookupTable.csv",
                                             mywd = "/data/notebooks/jupyterlab-jlconnectivity")})

sparkInputs

# start spark session 
sparkR.session(appname = "SparkR-Test",
               sparkHome = Sys.getenv("SPARK_HOME"),
               sparkConfig = list(spark.executor.instances = "3",
                                  spark.executor.cores = "4",
                                  spark.executor.memory = "8g",
                                  spark.kubernetes.executor.limit.cores = "4",
                                  spark.kubernetes.container.image = "nerc/sparkr-k8s:0.3.0"))

spark.lapply(sparkInputs, studyareabasics)

sparkR.session.stop()


# running functions for the whole of England

# system.time({england <- studyareabasics(regionfname = "England_AL4-AL4", 
#                                         regionspath = "/data/input-data/ToShare/Eng_shp", lcm = lcm,
#                                         resn = 25, align = 1000, hablabel = "heath",
#                                         habtext = "heathland", 
#                                         hlookuppath = "/data/input-data/ToShare/Habitat_lookupTable.csv")})

  