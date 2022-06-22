# SPATIAL USE NETWORKS

# This document contains the novel code of the manuscript titled
# 'Resource predictability modulates space use networks in an endangered scavenger species'
# by Cerecedo-Iglesias et al. submitted to Ecological Applications in 2022.

# Dealing with directory -------------------------------------------------------

# Creating working folders

dir <- getwd() # If required, set your directory with setwd() 

loc.data <- paste0(dir, "/DATA/")
loc.source <- paste0(dir, "/Functions/")

# Functions needed for the analysis are located in the 'Functions_SpatialNetworks.R' file
source(paste0(loc.source, 'Functions_SpatialNetworks.R'))

# Data not available (sensitive species) ---------------------------------------

# Load the cleaned data frame of locations of all individuals during the study period

load("all.data.RData")

# Calculating dynamic Brownian Bridge Models (dbbmm)----------------------------
# We calculate dbbmm by individual + year 

# Set the parameters

# specify the projection system in use
crs.proj <- CRS("+init=epsg:25831")

# Move::brownian.motion.variance.dyn parameters
ws <- 31
mrg <- 11
location.error= 30

# Move::brownian.bridge.dyn parameters
ext <- 0.5
r<-raster(xmn=-542632.6, xmx=740864.3, ymn=3842242, ymx=5093342, 
          crs=crs.proj, res=c(500,500))

# We have a gappy data because GPS are off during night.
# So we will ignore for all segments that have a larger time lag than 11.8 hours. 
# Therefore we consider all segments with time lag larger than 708 mins.
gappy_lag = 708 

dbmm_list <- dBBMM_calculation(all.data, r, ws, mrg, location.error, crs.proj, 
                                grid.ext, gappy_lag=gappy_lag, 
                                save_poly=FALSE)

# WARNING!!! We get a dbbmm by individual and year. As we need the information 
# at population level we calculate the average of all dBBMM layers.

# Calculating home ranges at population level ----------------------------------

# To convert dbbmm to rasterstack, we used the function developed by Virginia Morera-Pujol
# The function is available in the link: 
# https://github.com/VirginiaMorera/Useful-little-functions/blob/master/list_to_stack.R

dbbmm.list <- list_to_stack(dbbmm.list, 500, "+init=epsg:25831")
average_UD <- calc(dbbmm.list, fun = mean, na.rm = T)

# Extracting nodes at population level -----------------------------------------
# Nodes are defined as each isolated polygon of 50% countour home ranges

average_UD <- as(average_UD, 'DBBMM') 
poly50 <- raster2contour(average_UD, level=.50)
poly50 <- st_transform(st_as_sf(poly50, wkt = "geometry", crs=25831), crs=25831)
poly50<-st_make_valid(st_cast(poly50, 'POLYGON'))
poly50 <- st_cast(poly50, "POLYGON")


# Save separately all polygons of the home range in a list, hereafter nodes
nodes.list <- vector(mode='list')

for (j in 1:nrow(poly50)){
  poly_dis<-poly50[j, ]
  nodes.list[[j]] <- poly_dis
}

# Getting the nodes ------------------------------------------------------------

# First, calculate the centroids of nodes
nodes <- vector(mode = "list")
iter <- 0

for (i in 1:length(nodes.list)) {
  centroid <- nodes.list[[i]] %>%
    st_centroid() %>%
    st_geometry()
  
  iter <- iter +1
  nodes[[iter]] <- st_coordinates(centroid)
}

nodes <- as.data.frame(do.call(rbind, nodes))
colnames(nodes) <- c("long", "lat")
rownames(nodes) <- NULL

# Change nodes projection

nodes <- st_as_sf(nodes, coords = c("long", "lat"), crs = 25831)
nodes <- st_transform(nodes, crs = 4326)

# Getting trips ----------------------------------------------------------------

# First, annotation node by location

all.data <- st_as_sf(all.data, coords = c("long", "lat"), crs = 4326)

inter.all<-data.frame()

for(i in length(poly50)){
  
  print(paste("-------- node", i))
  
  s <- sf::st_transform(poly50[i, ], crs = 4326)
  
  inter <- sf::st_intersection(all.data, s)
  st_geometry(inter)<-NULL
  inter<-inter[1:4]
  inter$node <- paste("node", site[i])
  
  inter.all<- rbind(inter.all, inter)
}

load("all.data.RData") # Reload the clean dataframe

split.data <- merge(all.data, inter.all, by=c("id_device", "year", "DateTime", "id"), all.x=T)
split.data$node[is.na(split.data$node) == T]<-"none"

# Calculate trips between nodes ------------------------------------------------

trips <- trip_calculation(split.data, 20000)

# Trip 1 hour restriction

trips %>% group_by(track) %>% summarize( 
  ID            = unique(ID),
  DateTime_from = unique(DateTime_from),
  DateTime_to   = unique(DateTime_to),
  from          = unique(from),
  to            = unique(to),
  distance      = sum(distDiff, na.rm = T),
  duration      = sum(timeDiff, na.rm = T), # minutes
  duration_days = sum(timeDiff, na.rm = T)/(60*24), 
  angle_mean    = mean(angle, na.rm = T),
  speed_mean    = mean(speed, na.rm =T),
) %>% 
  filter(duration > 60) %>%
  droplevels() %>%
  arrange(ID, DateTime_from)-> trips # this contains all trips longer than 1 hours


# Here, we characterized the nodes (see Methods of main manuscript). We reload
# nodes data frame.

load('nodes.RData')

# Here, we show the function to print a spatial use network
# Must be a covariate called resource in nodes data frame

a <- network_graph(trips, nodes, plot = TRUE) 




