# Functions 

if(! 'move' %in% rownames(installed.packages())) install.packages('move')
if(! 'rgdal' %in% rownames(installed.packages())) install.packages('rgdal')
if(! 'sf' %in% rownames(installed.packages())) install.packages('sf')
if(! 'sp' %in% rownames(installed.packages())) install.packages('sp')
if(! 'igraph' %in% rownames(installed.packages())) install.packages('igraph')
if(! 'ggplot2' %in% rownames(installed.packages())) install.packages('ggplot2')
if(! 'dplyr' %in% rownames(installed.packages())) install.packages('dplyr')
if(! 'forecast' %in% rownames(installed.packages())) install.packages('forecast')
if(! 'marmap' %in% rownames(installed.packages())) install.packages('marmap')
if(! 'omxr' %in% rownames(installed.packages())) install.packages('omxr')
if(! 'raster' %in% rownames(installed.packages())) install.packages('raster')

library(move)
library(rgdal)
library(sf)
library(sp)
library(dplyr)
library(ggplot2)
library(igraph)
library(marmap)
library(forecast)
library(omxr)
library(raster)

dBBMM_calculation <- function(data, r, ws, mrg, location.error, crs.proj, 
                              grid.ext, gappy_lag=0, save_poly=FALSE) {
  #' Calculate dynamic Brownian Bridge Mododels (dBBMM)
  #'
  #' @description This function calculate the dBBMM from dataframe and 
  #' the kernel surface
  #' @param data dataframe with id, lat, long (WGS84) and Timestamp. 
  #' @param r raster See move::move help
  #' @param ws window size See move::move help
  #' @param mrg margin See move::move help
  #' @param location.error location.error See move::move help
  #' @param crs.proj crs projection of UDs and kernels
  #' @param ext See move::move help
  #' @param gappy_data defaults 0, if there is a gap in the data indicate the
  #' time lag when motion variance will be removed
  #' 
  #' @return two objects. First, the dataframe with kernel surface. Second,
  #' the dBBMM list
  
  # Iterations per individual and year: id = id_device + year
  who <- levels(as.factor(all.data$id))
  
  # dBBMM Loop per Individual --------------------------------------------------
  
  iter <- 0
  all.area.list <- vector(mode = "list", length = length(unique(who)))
  poly.list <- vector(mode = "list", length = length(unique(who)))
  dbbmm.list <- vector(mode = "list", length = length(unique(who)))
  
  for(indi in unique(who)){
    
    # pull out one individual
    data <- all.data[all.data$id == indi,]
    
    #iter <- iter + 1
    
    print(paste('----', indi, ': Creating move object'))
    
    # convert data into a move object ready for dBBMM calculation
    # Calculating trajectory:
    move <- move(x=data$long,
                 y=data$lat,
                 time = as.POSIXct(data$DateTime,format="%Y-%m-%d %H:%M:%s", 
                                   tz="UTC"),
                 proj=CRS("+proj=longlat +ellps=WGS84"),
                 animal=data$id, 
                 data=data,
                 removeDuplicatedTimestamps=TRUE)
    move <- spTransform(move, crs.proj)
    
    
    # Calculating dBBMM:
    ## calculate the dynamic brownian motion variance of the gappy track
    
    print(paste('----', indi, ': Calculating motion variance'))
    dbbv <- brownian.motion.variance.dyn(move, location.error=location.error, 
                                         window.size = ws, margin = mrg) 
    
    # if you have a gappy data
    
    dbbv@interest[timeLag(dbbv,"mins") > gappy_lag] <- FALSE 
    
    ## then we use the 'dBMvariance' object to calculate the dBBMM
    
    print(paste('----', indi, ': Calculating dBBMM'))
    dbbmm <- brownian.bridge.dyn(dbbv, raster = r, ext = ext, 
                                 location.error=30, verbose = FALSE)
  }
  
  dbbmm.list[[iter]] <- dbbmm
  
  return (dbbmm.list)
}


list_to_stack <- function(raster_list, new_res, dest_crs = CRS("+proj=longlat"), turn_0_to_NA = FALSE) {
  # Function developed by Virginia Morera-Pujol
  # The function is available in the link: 
  # https://github.com/VirginiaMorera/Useful-little-functions/blob/master/list_to_stack.R
  
  if(turn_0_to_NA == TRUE) {
    for(i in seq_along(raster_list)) {
      raster_list[[i]][raster_list[[i]] == 0] <- NA
      raster_list[[i]] <- trim(raster_list[[i]])}
  }
  
  # reproject rasters to dest_crs
  raster_list_proj <- lapply(raster_list, projectRaster, crs = dest_crs)
  
  # get the minimum common extent that encompasses all rasters in list 
  ext_list <- lapply(X = raster_list_proj, function(x) {as.matrix(x@extent)})
  matrix_extent <- matrix(unlist(ext_list), ncol=length(ext_list))
  rownames(matrix_extent)<-c("xmin", "ymin", "xmax", "ymax")
  best_extent <- extent(min(matrix_extent[1,]), max(matrix_extent[3,]), 
                        min(matrix_extent[2,]), max(matrix_extent[4,]))
  
  # resample all rasters in list to best extent, new crs and new resolution
  res_list2 <- lapply(raster_list_proj, resample, 
                      raster(ext = best_extent, crs = dest_crs, resolution = new_res))
  
  # turn to stack
  res_stack <- stack(res_list2)
  
  return(res_stack)
}



network_graph <- function(data, nodes, node_feature = NULL, plot = TRUE){
  
  #' Plot a spatial use network
  #'
  #' @description This function plot a spatial use network
  #' @param data data frame with id, from node and to node. It dataframe must to be all
  #' trips between nodes. Names of column must to be 'id', 'from' and 'to'
  #' @param nodes dataframe with nodes names and their features. The firs column
  #' must to be the node id and must to call 'nodes'
  #' @param node_feature covariate which characterized the nodes.
  #' @param plot logical. If TRUE, plot he graph
  #' 
  #' @return a plot
  
  # Iterations per individual and year: id = id_device + year
  
  # Getting from and to nodes 
  from <- data$from
  to <- data$to
  directions <- data.frame(from=from, to=to)
  
  P <-as.data.frame.matrix(table(directions))
  P <- gather_matrix(as.matrix(P), value_name="weight")
  colnames(P) <- c("from", "to", "weight") 
  P <- subset(P, weight > 0 )
  
  net <- graph_from_data_frame(d=P, vertices=nodes, directed=T)
  local.metrics.nobreeders <- data.frame(edges = gsize(net), 
                                         ver = length(V(net)[degree(net)>0]),
                                         bet = betweenness(net), deg = degree(net), 
                                         degin = degree(net, mode = "in"), 
                                         degout = degree(net, mode = "out"))
  local.metrics.nobreeders$nodes <- rownames(local.metrics.nobreeders)
  rownames(local.metrics.nobreeders)<-NULL
  
  N<-merge(nodes, local.metrics.nobreeders, by = ("nodes"))
  
  edges_for_plot <- P %>% 
    inner_join(N %>% dplyr::select(nodes, long, lat), by = c('from' = 'nodes')) %>%
    rename(xfrom = long, yfrom = lat) %>%
    inner_join(N %>% dplyr::select(nodes, long, lat), by = c('to' = 'nodes')) %>%
    rename(xto = long, yto = lat)
  
  edges_for_plot <- subset(edges_for_plot, xfrom != xto)
  
  land <- getNOAA.bathy(lon1 = -3, lon2 = 5,
                        lat1 = 40, lat2 = 44, resolution = 1) 
  l <- forecast::autoplot(land, geom=c("r"), show.legend = F) + 
    scale_fill_etopo() 
  
  q <- l + 
    geom_curve(aes(x = xfrom, y = yfrom, xend = xto, yend = yto,    
                   size = weight/10),
               data = edges_for_plot, curvature = 0.33,
               alpha = 0.5) +
    scale_size_continuous(guide = FALSE, range = c(1, 8)) +
    ylab("lat") + xlab("long") 
  
    a <- q + geom_point(data=N, aes(x = long, y = lat, size = bet, 
                                    color = resource), alpha = 0.8) +
      scale_color_brewer(palette="Dark2") +
      theme_bw() + 
      ggtitle("Spatial use network") + 
      theme(axis.text = element_text(size = 9),
            axis.title.y = element_text(size=9),
            plot.title = element_text(size=10))
  
  if (plot == TRUE) {
    plot(a)
  }
  return(a)
}



random_remove <- function(nodes, freq_mat, number_removal, thershold, iterations) {
  
  #' Calculate the global measures and nodes measures of random removal simulations
  #'
  #' @description This function simulates the removal of nodes randomly
  #' @param nodes data frame with nodes names and their features.
  #' The first column must  be the node id and must call 'nodes'
  #' @param freq_data data frame with 'from' node and 'to' node. 
  #' This data frame must incorporate all trips between nodes. Names of columns must be 'from' and 'to'.
  #' @param number_removal Number of nodes to be removed
  #' @param threshold limit of nodes to remove. Values must be between 0-1.
  #' @param iterations number of iteration
  #' 
  #' @return a list of data frames with global measures and node measures in each 
  #' iteration of modelling.

  iter <- 0
  global.measures <- vector(mode = "list")
  nodes.measures <- vector(mode = "list")
  number_removal <- seq(1, number_removal, 1)
  
  for(i in number_removal){ 
    p <- number_removal[i]/nrow(nodes) 
    
    if(p < thershold){
      
      print(paste("------", i, "-----", p))
      
      for(j in 1:iterations){ 
        
        rdm <- sample(nodes$nodes, size = i) 
        
        net <- graph_from_data_frame(d=freq_mat, vertices=nodes, directed=T)
        for (name in rdm) {
          net <- delete_vertices(net, name)
        }
        
        
        global.metrics.pop <- data.frame(rem = i, rep = j, den = edge_density(net),
                                         avl = mean_distance(net), trans = transitivity(net), 
                                         diam = diameter(net))
        
        local.metrics.pop <- data.frame(rem = i, rep = j, bet = betweenness(net))
        
        iter <- iter + 1
        
        global.measures[[iter]] <- global.metrics.pop
        nodes.measures[[iter]] <- local.metrics.pop
      }
    } else {print('finish')
      break}
  }
  
  global.measures.random <- as.data.frame(do.call(rbind, global.measures))
  global.measures.random$removal <- rep("random", nrow(global.measures.random))
  nodes.measures.random <- as.data.frame(do.call(rbind, nodes.measures))
  nodes.measures.random$removal <- rep("random",  nrow(nodes.measures.random))
  
  return(list(global.measures.random, nodes.measures.random))
  
}



target_remove <- function(nodes, target, freq_mat, number_removal, threshold, iterations) {
  
  #' Calculate the global measures and nodes measures of targeted removal simulations
  #'
  #' @description This function simulates the removal of nodes of a spatial use 
  #' network based on a specific feature of nodes
  #' @param nodes data frame with nodes names and their features. 
  #' The first column must  be the node id and must call 'nodes'
  #' @param target subset of nodes data frame focus to remove
  #' @param freq_data data frame with 'from' node and 'to' node. 
  #' This data frame must incorporate all trips between nodes. Names of columns must be 'from' and 'to'.
  #' @param number_removal Number of nodes to be removed
  #' @param threshold limit of nodes to remove. Values must be between 0-1.
  #' @param iterations number of iteration
  #' 
  #' @return a list of data frames with global measures and node measures in each 
  #' iteration of modelling.
  
  iter <- 0
  global <- vector(mode = "list")
  global.measures <- vector(mode = "list")
  nodes.measures <- vector(mode = "list")
  number_removal <- seq(1, number_removal, 1)
  
  for(i in number_removal){ 
    p <- number_removal[i]/nrow(nodes) 
    
    if(p < threshold){
      
      print(paste("------", i, "------", p))
      
      for(j in 1:iterations){ 
        
        rdm <- sample(target$nodes, size = i) 
        
        net <- graph_from_data_frame(d=freq_mat, vertices=nodes, directed=T)
        for (name in rdm) {
          net <- delete_vertices(net, name)
        }
        
        global.metrics.pop <- data.frame(rem = i, rep = j, den = edge_density(net),
                                         avl = mean_distance(net), trans = transitivity(net), 
                                         diam = diameter(net))
        
        local.metrics.pop <- data.frame(rem = i, rep = j, bet = betweenness(net))
        
        iter <- iter + 1
        
        global.measures[[iter]] <- global.metrics.pop
        nodes.measures[[iter]] <- local.metrics.pop
      }
    } else {print('finish') 
      break}
  }
  
  global.measures.target <- as.data.frame(do.call(rbind, global.measures))
  global.measures.target$removal <- rep("target", nrow(global.measures.target))
  nodes.measures.target <- as.data.frame(do.call(rbind, nodes.measures))
  nodes.measures.target$removal <- rep("target",  nrow(nodes.measures.target))
  
  return(list(global.measures.target, nodes.measures.target))
  
}


