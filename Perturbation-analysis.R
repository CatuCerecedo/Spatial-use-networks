# PERTURBATION ANALYSIS

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

# Load data --------------------------------------------------------------------

load('nodes.RData')
load('trips.RData')

# Perturbation analysis --------------------------------------------------------

# Random remove
global.measures.random <- random_remove(nodes, 
                                        trips, 41, 0.5, 1000)
nodes.measures.random <- global.measures.random[[2]]

# Target 
# We remove nodes by resources availability feature

### Landfill removal 

target <- nodes %>% filter (resource == "Landfill") 

global.measures.target <- target_remove(nodes, target, trips, 7, 0.5, 1000)
nodes.measures.target <- global.measures.target[[2]] ###### De donde sale esto!!!! Comrpobar función

# Save simulations measures

nodes.measures<- rbind(nodes.measures.target, nodes.measures.random)

data <- nodes.measures %>% filter(rem < 8) %>% 
  transform(rem = as.factor(rem)) 

a<-ggline(data, x="rem", y="bet", add="mean_ci", size = 1, point.size = 1,
          shape ="removal", linetype = "removal", position=position_dodge(0.1))  +
  stat_compare_means(aes(group = removal), hide.ns = TRUE, method = "t.test",
                     paired = T, label = "p.signif", symnum.args = psymbols,
                     size = 5, label.y = 83) +
  xlab("") + ylab("Non-breeders\n\nBetweenness") +
  theme_bw() + ggtitle("a) Landfill") +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size=9),
        plot.title = element_text(size=10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10)) 


### Vultures restaurant removal 

target <- nodes %>% filter (resource == "Vulture restaurant") 

global.measures.target <- target_remove(nodes, target, trips, 3, 0.5, 1000)
nodes.measures.target <- global.measures.target[[2]]

# Save simulations measures

nodes.measures<- rbind(nodes.measures.target, nodes.measures.random)

data <- nodes.measures %>% filter(rem < 4) %>% 
  transform(rem = as.factor(rem))


b<-ggline(data, x="rem", y="bet", add="mean_ci", size = 0.6, point.size=1,
          shape ="removal", linetype = "removal", position=position_dodge(0.1))  +
  stat_compare_means(aes(group = removal), hide.ns = TRUE, method = "t.test",
                     paired = T, label = "p.signif", symnum.args = psymbols,
                     size = 5, label.y = 84.7) +
  xlab("") + ylab("") +
  theme_bw() + ggtitle("b) Vultures restaurant") +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size=9),
        plot.title = element_text(size=10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10)) 

### Intensive farm removal 

target <- nodes %>% filter (resource == "Intensive") 

global.measures.target <- target_remove(nodes_nobreeders, target, P.nobreeder,5, 0.5, 1000)
nodes.measures.target <- global.measures.target[[2]]

# Save simulations measures

nodes.measures<- rbind(nodes.measures.target, nodes.measures.random)

data <- nodes.measures %>% filter(rem < 5) %>% 
  transform(rem = as.factor(rem))

c<-ggline(data, x="rem", y="bet", add="mean_ci", size = 0.6, point.size=1,
          shape ="removal", linetype = "removal", position=position_dodge(0.1))  +
  stat_compare_means(aes(group = removal), hide.ns = TRUE, method = "t.test",
                     paired = T, label = "p.signif", symnum.args = psymbols,
                     size = 5, label.y = 83) +
  xlab("") + ylab("") +
  theme_bw() + ggtitle("c) Intensive farm") +
  theme(axis.text = element_text(size = 9),
        axis.title.y = element_text(size=9),
        plot.title = element_text(size=10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10)) 

### Extensive farm removal 

target <- nodes %>% filter (resource == "Extensive") 

global.measures.target <- target_remove(nodes, target, trips, 10, 0.5, 1000)
nodes.measures.target <- global.measures.target[[2]]

# Save simulations measures

nodes.measures<- rbind(nodes.measures.target, nodes.measures.random)

data <- nodes.measures %>% filter(rem < 11) %>% 
  transform(rem = as.factor(rem))

d<-ggline(data, x="rem", y="bet", add="mean_ci", size = 0.6, point.size=1,
          shape ="removal", linetype = "removal", position=position_dodge(0.1))  +
  stat_compare_means(aes(group = removal), hide.ns = TRUE, method = "t.test",
                     paired = T, label = "p.signif", symnum.args = psymbols,
                     size = 5, label.y = 83) +
  xlab("") + ylab("") +
  theme_bw() + ggtitle("d) Extensive farm") +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size=9),
        plot.title = element_text(size=10),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10)) 

# All plot together
ggarrange(a,b,c,d, nrow = 1, ncol = 4, common.legend = T, legend = "right")
