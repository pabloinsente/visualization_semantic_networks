# Author: Pablo Caceres
# Created: 03/29/2019
# Updated: 04/29/2019
# R-version: 3.5.1 

########################
# Next data preparation
########################

###############################################################################################
# This section:
# - Reads in the data frame
# - Creates a unique labels' dataframe to feed into the graph object as the nodes/vertices list
# - Creates a link's dataframe: 
#                               - "Center" is taken as the "From" node
#                               - "Answer" is taken as the "To" node
#                               - All nodes are directed
################################################################################################

# load data
d_raw <- read.csv(file="data/animal_random_merged.csv", header=TRUE, sep=",") 
# create node list of unique labels
nodes <- data.frame(unique(d_raw$Center)) 
# change colname Label
nodes <- setNames(nodes, c("Label")) 
# vectors to pick
myvars <- c("Center", "Answer") 
# pick "Center" and "Answer" vectors 
links <- d_raw[myvars]  
# add vector "Type" of link 
links$Type <- rep(c("Directed"), times = length(links$Center)) 
# change col names for graph object
links <- setNames(links, c("From", "To", "Type")) 

###########################
# Creating a Network Object
###########################

library(igraph)

# we need to pass "link" as edges to "d" and "nodes" as vertices to "Vertices"
# we also set "directed" to true ("T")
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) # new igraph-network object
E(net) # edges
V(net) # vertices

#####################
# Five-number summary
#####################
library(intergraph)
library(statnet)

# We first need to coerce the igraph object to a statnet network object
net.n <- asNetwork(net)

# Size: number of nodes
network.size(net.n)
## Size: 217 ##

# Density: proportion of observed ties of the maximum number of possible ties.
gden(net.n)
## Density: 0.21178955453149 ##

# Components: is a subgroup in which all actors are connected, directly or indirectly. 

#Find strong components
comp = components(net.n,connected="strong")
# check help(component.dist) for details about how components are identified
## Components: 2 ##

# Diameter (geodesic distance): is the longest of the shortest paths across all pairs of nodes
# Is measure of compactness
lgc <- component.largest(net.n,result="graph")
gd <- geodist(lgc)
max(gd$gdist)
## Diameter: 3 ##

# Clustering coefficient: is the presence of clustering, or the tendency to formed closed triangles

#we first need to check if the network is multiplex 
# (multiple directed edges between a pir of nodes)
is.multiplex(net.n)

# It is, so let's remove loops and multiple edges to simplify the net
net.s <- simplify(net, remove.multiple = T, remove.loops = T)

# Let's coerce the igraph object to a statnet network object again
net.n.s <- asNetwork(net.s)

# Now use gtrans to compute diamater
gtrans(net.n.s,mode="graph")
## Clustering coefficient: 0.198555000965171 ##

###############
# Basic Summary
###############
## Size:       217  ##
## Density:    0.21 ##
## Components: 2    ##
## Diameter:   3    ##
## Clustering: 0.19 ##

########################
# Basic Network plotting 
########################
library(scales)
library(beepr)

# First raw plot
plot(net)

# Compute node degrees (#links) and use that to set node size:
deg <- degree(net.n.s)
deg <- rescale(deg, to = c(5, 15)) # resize node size 
V(net.s)$size <- deg
plot(net.s, edge.arrow.size=.1) # reduce arrow size

# Let's remove labels from the graph to see how it looks
lk = layout.kamada.kawai(net.s) # create a tree like layout
plot(net.s, edge.arrow.size=.1,vertex.label=NA, layout = lk)
# Really clumpy!

#####################
# Community detection
#####################

# cluster edge betweenness
# The idea of the edge betweenness based community structure detection 
# is that it is likely that edges connecting separate modules have high 
# edge betweenness as all the shortest paths from one module to another 
# must traverse through them.

# This may take a minute to run...
clp_e <- cluster_edge_betweenness(net.s)
beep(sound = 8) # sound when finish

# cluster random walk
# This function tries to find densely connected subgraphs, also called 
# communities in a graph via random walks. The idea is that short random
# walks tend to stay in the same community.

# Note: this works only for undirected graphs, hence, it will ignore 
# edge direction

clp_r <- cluster_walktrap(net.s)
beep(sound = 2) # sund when finish

# The cluster function returns a community object
class(clp_e)
class(clp_r)

# number of communities detected by each algorithm
length(clp_e) 
length(clp_r) 

# here we can check the community at which each node was assigned 
head(clp_e$membership) 
head(clp_r$membership) 

# how modular the graph partitioning is
# High modularity for a partitioning reflects dense connections within 
# communities and sparse connections across communities.

modularity(clp_e)
modularity(clp_r)

###############
# Basic Summary
###############
# Edge betweenness ##
# Communities: 24
# Modularity: 0.00067
#####################
# Random trap (walk)#
# Communities: 4
# Modularity: 0.07381
#####################

# IMPORTANT: There is massive difference between detecting communties
# by a random walk that assumes undirected edges, and edge_betweenness
# that works with directed edges

#################
# Better plotting
#################

# create a tree like layout
lk = layout.kamada.kawai(net.s) 

# Let's create a color pallete for edge_betweenness communities
prettyColors <- rainbow(length(clp_e), alpha = 0.7)
communityColors_e <- prettyColors[membership(clp_e)]

# Let's create a color pallete for random_walk_trap communities
prettyColors <- rainbow(length(clp_r), alpha = 0.7)
communityColors_r <- prettyColors[membership(clp_r)]

### GRAPH 1 ##
# Let's plot the edge_betweenness result with lk layout
# Note: since there are communities of "1" or very few nodes, the layout.fruchterman.reingold
# won't work well

set.seed(28) # set seed for reproducibility

plot(clp_e, net.s, edge.width = 0.1, mark.border=NA, vertex.label.cex=0.6, 
     vertex.label.color="blue", vertex.size = 10, edge.arrow.size=0.1, 
     edge.color = "Grey", col = communityColors_e, layout= lk, asp =0)

### GRAPH 2 ##
# Let's plot the random_walk_trap result with lk layout
# Note: since there are communities of "1" or very few nodes, the layout.fruchterman.reingold
# won't work well

set.seed(28) # set seed for reproducibility

plot(clp_r, net.s, edge.width = 0.1, mark.border=NA, vertex.label.cex=0.6, 
     vertex.label.color="blue", vertex.size = 10, edge.arrow.size=0.1, 
     edge.color = "Grey", col = communityColors_r, layout= lk, asp =0)

## Plot out GRAPH 1 ##
set.seed(28) # set seed for reproducibility

tiff("plot_next_animal_eb.tiff", width =5, height = 5, units = 'in', res = 300, compression = 'none')

plot(clp_e, net.s, edge.width = 0.1, mark.border=NA, vertex.label.cex=0.5, 
     vertex.label.color="blue", vertex.size = 10, edge.arrow.size=0.1, 
     edge.color = "Grey", col = communityColors_e, layout= lk, asp =0)

dev.off()

## Plot out GRAPH 2 ##
set.seed(28) # set seed for reproducibility

tiff("plot_next_animal_rw.tiff", width =5, height = 5, units = 'in', res = 300, compression = 'none')

plot(clp_r, net.s, edge.width = 0.1, mark.border=NA, vertex.label.cex=0.5, 
     vertex.label.color="blue", vertex.size = 10, edge.arrow.size=0.1, 
     edge.color = "Grey", col = communityColors_r, layout= lk, asp =0)

dev.off()