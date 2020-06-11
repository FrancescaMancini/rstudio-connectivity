############################################
## Data manipulation
############################################

library(rgdal)

data_path <- "/data/input-data"


## extract single counties polygons from all counties shapefile

counties <- readOGR("/data/input-data/ToShare/Eng_counties", 
                    "Eng_Boundary-line-ceremonial-counties_region")

#West Yorkshire

w_yorkshire <- subset(counties, NAME == "West Yorkshire")

plot(w_yorkshire)

writeOGR(w_yorkshire, "/data/derived", "w_yorkshire", driver="ESRI Shapefile")

# West Sussex

w_sussex <- subset(counties, NAME == "West Sussex")

plot(w_sussex)

writeOGR(w_sussex, "/data/derived", "w_sussex", driver="ESRI Shapefile")
