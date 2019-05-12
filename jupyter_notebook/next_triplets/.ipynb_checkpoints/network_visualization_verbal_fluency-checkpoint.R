# Author: Pablo Caceres
# Created: 03/29/2019
# Updated: 05/12/2019
# R-version: 3.6.0

# load libraries
library(igraph)
library(scales)
library(lmSupport)
library(beepr)
library(RColorBrewer)
library(qgraph)

#################################
# Verbal Fluency data preparation
#################################

# Read in "animal" data
nodes <- read.csv("data/animal/animals_min5_directed_Nodes.csv", header=T, as.is=T)
links <- read.csv("data/animal/c05_st_red_anim_min5_dir_edges.csv", header=T, as.is=T)

# Examine the data
head(nodes)
head(links)

###########################
# Creating a Network Object
###########################

# Converting the data to an igraph object
# Arguments:
# - d = edges (links)
# - vertices = nodes
# - directed = True (if directed graph)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
is_weighted(net) # check if is weigthed

# Examine the resulting object
class(net)
net 

# We can access the nodes, edges, and their attributes:
E(net) # edges
V(net) # vertices
print("Attributes edges")
head(E(net)$type) # attributes edges
print("Attributes vertex")
head(V(net)$name) # attributes vertex 

# To access specific nodes and edges by attribute
V(net)[name=="aardvark"]
E(net)[type=="Directed"]

# To extract an edge list or a matrix back from the igraph networks:
head(as_edgelist(net, names=T))
head(as_adjacency_matrix(net, attr="weight"))

# To extract data frames describing nodes and edges:
head(as_data_frame(net, what="edges"))
head(as_data_frame(net, what="vertices"))

#####################
# Five-number summary
#####################

library(intergraph)
library(statnet)

# We first need to coerce the igraph object to a statnet network object
net.n <- asNetwork(net)

# Size: number of nodes
network.size(net.n)
## Size: 295 ##

# Density: proportion of observed ties of the maximum number of possible ties.
gden(net.n)
## Density: 0.11293670010377 ##

# Diameter (geodesic distance): is the longest of the shortest paths across all pairs of nodes
# Is measure of compactness
lgc <- component.largest(net.n,result="graph")
gd <- geodist(lgc)
max(gd$gdist)
## Diameter: 4 ##

# Components: is a subgroup in which all actors are connected, directly or indirectly. 

#Find strong components
comp = components(net.n,connected="strong")
# check help(component.dist) for details about how components are identified
comp
## Components: 1 ##

# Clustering coefficient: is the presence of clustering, or the tendency to formed closed triangles

#we first need to check if the network is multiplex 
# (multiple directed edges between a pir of nodes)

is.multiplex(net.n) # FALSE

# It is not, so we can procede

# Now use gtrans to compute diamater
gtrans(net.n,mode="graph")
## Clustering coefficient: 0.441780172888112 ##

###############
# Basic Summary
###############
## Size:       295  ##
## Density:    0.11 ##
## Components: 1    ##
## Diameter:   4    ##
## Clustering: 0.44 ##

##########
# Plotting
##########

# First raw plot
plot(net)

# Simplified network removing loops and multiple edges between two nodes
net.s <- simplify(net, remove.multiple = T, remove.loops = T, edge.attr.comb=c(weight="sum", type="ignore"))
plot(net.s)

# Compute node degrees (#links) and use that to set node size:
deg <- degree(net.s, mode="all")
deg <- rescale(deg, to = c(4, 20)) # resize node zies between 4-20
V(net.s)$size <- deg
plot(net.s, edge.arrow.size=.1) # reduce arrow size

# Let's remove labels from the graph to see how it looks
plot(net.s, edge.arrow.size=.1,vertex.label=NA)

# Let's reduce arrow size, put labels back, and reuce label size 
plot(net.s, edge.arrow.size=.1, vertex.color="orange", 
     vertex.frame.color="#555555",vertex.label=V(net)$name, 
     vertex.label.color="black", vertex.label.cex=.5)

# Let's try removing the vertex shape (circles), changing labels to blue
plot(net.s, vertex.shape="none", vertex.label=V(net)$name, 
     vertex.label.font=2, vertex.label.color="blue",
     vertex.label.cex=.5, edge.arrow.size=.1,
     edge.color="gray85")

# Let's try differnet layouts by using the layout plot function:
# ?igraph::layout_ To check all layouts

# randomly placed vertices
plot(net.s, vertex.label=V(net)$name, vertex.label.font=2, 
     vertex.label.color="blue", vertex.label.cex=.5, 
     edge.arrow.size=.1, edge.color="gray85",layout=layout_randomly)

# 3D sphere layout
plot(net.s, vertex.label=V(net)$name, vertex.label.font=2, 
     vertex.label.color="blue",vertex.label.cex=.5, 
     edge.arrow.size=.1, edge.color="gray85", layout=layout_on_sphere)

#lgl layout works better with big-dense networks
plot(net.s, vertex.label=V(net)$name, vertex.label.font=2, 
     vertex.label.color="blue",vertex.label.cex=.5, 
     edge.arrow.size=.1, edge.color="gray85", layout=layout_with_lgl)

#####################################
# Highlighting aspects of the network
#####################################

# Lte's look at the weights' distribution
hist(links$weight, breaks=seq(1, max(links$weight)))
varDescribe(links$weight) # describe weights
# Problem: we have way too many links and super skeweed

# To help to visualize the network, we are going to trim all the edges with N# of connections below the mean
cut.off <- mean(links$weight) 
net.sp <- delete_edges(net.s, E(net.s)[weight<cut.off])
plot(net.sp, vertex.label.color="blue", vertex.label.cex=0.5, 
     edge.arrow.size=0.1, edge.color="gray85",layout=layout_with_kk) 


### Community detection algorithms 
# We are going to try to make the network map more useful by 
# showing the communities within it.
# Community detection (by optimizing modularity over partitions) may take a 
# long with hundres of edgeds and  vertices  
# A few common options: 
# - **cluster_optimal**: This function calculates the optimal community structure of a graph, by maximizing the modularity measure over all possible partitions.    
# - **cluster_fast_greedy**: This function tries to find dense subgraph, also called communities in graphs via directly optimizing a modularity score.   
# - **cluster_walktrap**: This function tries to find densely connected subgraphs, also called communities in a graph via random walks. The idea is that short random walks tend to stay in the same community.   
#
# *cluster_optimal* may take forever to run with a large graph, and *cluster_fast_greedy* only works for undirected graphs, so so we are going to try with ***cluster_walktrap*** first.
# More info here https://stackoverflow.com/questions/9471906/what-are-the-differences-between-community- detection-algorithms-in-igraph

clp <- cluster_walktrap(net.s)
beep(sound = 8) # this will sound when its ready
# Super fast!

# The cluster_walktrap returns a community object
class(clp)
head(clp$membership) #here we can check the community at which each node was assigned 
print("number of communities: ")
length(clp) # number of communities
print("Graph's modularity: ")
modularity(clp) # how modular the graph partitioning is
#High modularity for a partitioning reflects dense connections within communities and sparse connections across communities.

# Let's look at the communities with a dendogram 
dendPlot(clp, mode="hclust")

# We can plot using the results obtained from the algorithm
# Igraph has built-in functions to plot this
#plot(clp, net, layout=layout_with_kk)
# But warning: PANDEMONIUM!
plot(clp, net.s, vertex.label.cex=0.5, edge.arrow.size=0.1,  
     edge.color="gray85",layout=layout_with_kk) 

# Lets create a color pallet for the 22 communities
prettyColors <- rainbow(22, alpha = 0.8)
communityColors <- prettyColors[membership(clp)]
plot(x=1:22, y=1:22, pch=19, cex=5, col=unique(communityColors))

# Now we can  plot the communities without relying on their built-in plot:
plot(clp, net.s, vertex.color=communityColors,vertex.label.cex=0.5,
     edge.arrow.size=0.1, edge.color="gray85",layout=layout_with_fr) 

# Function to group communities in the graph 
# Use:
# - put your community object in the community slot
# - put your graph in the network slot
# Mpre info: https://stackoverflow.com/questions/16390221/how-to-make-grouped-layout-in-igraph
edge.weights <- function(community, network, weight.within = 100, weight.between = 1) {
bridges <- crossing(communities = community, graph = network)
weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
return(weights) 
}

# Then transfer the weights to the weight slot in the nodes:
net.graph <- net.s
E(net.graph)$weight <- edge.weights(clp, net.graph)

# Plot with no labels, clustered by community, and nodes scaled by weight 
set.seed(28) # set seed for reproducibility
plot(x = clp, y = net.graph, edge.width = 0.1, mark.groups = NULL, ylim= c(-1, 1), xlim=c(-1,1), 
     vertex.label = NA, edge.arrow.size=0.1, edge.color = "Grey", col = communityColors, 
     layout = layout_with_fr, main = "Animal's association network")

# Plot with no labels, clustered by community, and nodes at fix size
set.seed(28) # set seed for reproducibility
plot(x = clp, y = net.graph, edge.width = 0.1,  mark.border=NA,
     vertex.label = NA, edge.arrow.size=0.1, vertex.size = 5, edge.color = "Grey",
     col = communityColors, layout = layout_with_fr,  main = "Animal's association network")

# Plot with labels, no shape,  clustered by community
set.seed(28) # set seed for reproducibility
plot(x = clp, y = net.graph, edge.width = 0.1, mark.border=NA, vertex.shape="none", vertex.label = V(net.graph)$name, 
     edge.arrow.size=0.1, vertex.label.font=1, edge.color = "Grey", vertex.label.color="black",
     vertex.label.cex=.6, col = communityColors, layout = layout_with_fr,  main = "Animal's association network")

# Options to plot out the graph in high resolution: https://www.r-bloggers.com/high-resolution-figures-in-r/

# Export as high resolution tiff

# Plot with no labels, clustered by community, and nodes scaled by weight 
# Export as high resolution tiff

set.seed(28) # set seed for reproducibility
tiff("plot_weighted_nodes_white_b.tiff", width =5, height = 5, units = 'in', res = 300, compression = 'none')
plot(x = clp, y = net.graph, edge.width = 0.1,  mark.border=NA,
     vertex.label = NA, edge.arrow.size=0.1, edge.color = "Grey",
     col = communityColors, layout = layout_with_fr)
dev.off()

# Note about plot layout
# By default, the coordinates of the plots are rescaled to the [-1,1] interval for both x and y. 
# We can change that with the parameter rescale=FALSE and rescale the plot manually by multiplying 
# the coordinates by a scalar. You can use norm_coords to normalize the plot with the boundaries we want.

l <- layout_with_fr(net.graph)

l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

par(mfrow=c(2,2), mar=c(0,0,0,0))

set.seed(28) # set seed for reproducibility
plot(x = clp, y = net.graph, edge.width = 0.1,  mark.border=NA,
     vertex.label = NA, edge.arrow.size=0.1, vertex.size = 5, edge.color = "Grey",
     col = communityColors, rescale=F, layout = l*0.4)
plot(x = clp, y = net.graph, edge.width = 0.1,  mark.border=NA,
     vertex.label = NA, edge.arrow.size=0.1, vertex.size = 5, edge.color = "Grey",
     col = communityColors, rescale=F, layout = l*0.6)
plot(x = clp, y = net.graph, edge.width = 0.1,  mark.border=NA,
     vertex.label = NA, edge.arrow.size=0.1, vertex.size = 5, edge.color = "Grey",
     col = communityColors, rescale=F, layout = l*0.8)
plot(x = clp, y = net.graph, edge.width = 0.1,  mark.border=NA,
     vertex.label = NA, edge.arrow.size=0.1, vertex.size = 5, edge.color = "Grey",
     col = communityColors, rescale=F, layout = l*1.0)

# Let's extract the first animal from each cluster to use it as label for the group
names = c()
for (i in 1:22) {
    names[[i]] <- paste(clp[[i]][1], ", n = ",sizes(clp)[[i]][1], sep = "")
    # append to list <- first animal in the cluster + number of animals in the cluster
}
names
class(names)

# GRAPH 1
# Let's add leyends to the groups by picking the first animal on each cluster
# Plot with no labels, clustered by community, and nodes at fix size

set.seed(28) # set seed for reproducibility
par(bg="black") # set background to black

plot(x = clp, y = net.graph, edge.width = 0.1,  mark.border=NA, 
     vertex.label = NA, edge.arrow.size=0.1, vertex.size = 5, edge.color = "Grey",
     col = communityColors, layout=layout_with_fr)

title(main = list("Animal's association network", cex = 1, col = "white", font = 2))

legend(x =-1.6 , y = 0.5, legend=levels(factor(names, levels=names)), 
       pt.bg=unique(communityColors), bty = "n", pch=20 , pt.cex = 1.0, cex = 0.6, 
       text.col=unique(communityColors), horiz = FALSE, inset = c(0.1, 0.1), text.font =2)

# Let's plot out graph 1
tiff("plot_unweighted_nodes_black_b.tiff", width =6, height = 5, units = 'in', res = 300, compression = 'none')

set.seed(28) # set seed for reproducibility
par(bg="black") # set background to black

plot(x = clp, y = net.graph, edge.width = 0.1,  mark.border=NA, 
     vertex.label = NA, edge.arrow.size=0.1, vertex.size = 5, edge.color = "Grey",
     col = communityColors, layout=layout_with_fr)

title(main = list("Animal's association network", cex = 1, col = "white", font = 2))

legend(x =-1.8 , y = 0.5, legend=levels(factor(names, levels=names)),
       col = unique(communityColors), bty = "n", pch=20 , pt.cex = 1.0, cex = 0.6, 
       text.col=unique(communityColors), horiz = FALSE, inset = c(0.1, 0.1), text.font =2)

dev.off()

# GRAPH 2
# Let's add leyends to the groups by picking the first animal on each cluster
# Plot with no labels, clustered by community, and nodes scaled by weight (n# of links)
set.seed(28) # set seed for reproducibility
par(bg="black") # set background to black
plot(x = clp, y = net.graph, edge.width = 0.1,  mark.border=NA, 
     vertex.label = NA, edge.arrow.size=0.2, edge.color = "Grey",
     col = communityColors, layout=layout_with_fr)

title(main = list("Animal's association network", cex = 1, col = "white", font = 2))

legend(x =-1.6 , y = 0.5, legend=levels(factor(names, levels=names)),
       col = unique(communityColors), bty = "n", pch=20 , pt.cex = 1.0, cex = 0.6, 
       text.col=unique(communityColors), horiz = FALSE, inset = c(0.1, 0.1), text.font =2)

# Let's plot out graph 2
tiff("plot_weighted_nodes_black_b.tiff", width =6, height = 5, units = 'in', res = 300, compression = 'none')

set.seed(28) # set seed for reproducibility
par(bg="black") # set background to black

plot(x = clp, y = net.graph, edge.width = 0.1,  mark.border=NA, 
     vertex.label = NA, edge.arrow.size=0.1, edge.color = "Grey",
     col = communityColors, layout=layout_with_fr)

title(main = list("Animal's association network", cex = 1, col = "white", font = 2))

legend(x =-1.8 , y = 0.5, legend=levels(factor(names, levels=names)),
       col = unique(communityColors), bty = "n", pch=20 , pt.cex = 1.0, cex = 0.6, 
       text.col=unique(communityColors), horiz = FALSE, inset = c(0.1, 0.1), text.font =2)
dev.off()

# GRAPH 3
# Let's add leyends to the groups by picking the first animal on each cluster
# Plot with labels, clustered by community, and nodes at fix size
set.seed(28) # set seed for reproducibility
par(bg="black") # set background to black

plot(x = clp, y = net.graph, edge.width = 0.1,  mark.border=NA,  vertex.label.cex=0.5, 
      vertex.label.color="white", edge.arrow.size=0.1, vertex.size = 6, edge.color = "Grey",
     col = communityColors, layout=layout_with_fr)

title(main = list("Animal's association network", cex = 1,col = "white", font = 2))

legend(x =-1.8 , y = 0.5, legend=levels(factor(names, levels=names)),
       col = unique(communityColors), bty = "n", pch=20 , pt.cex = 1.0, cex = 0.6, 
       text.col=unique(communityColors), horiz = FALSE, inset = c(0.1, 0.1), text.font =2)

# Let's plot out GRAPH 3
tiff("plot_unweighted_nodes__black_b_leg.tiff", width =6, height = 5, units = 'in', res = 300, compression = 'none')

set.seed(28) # set seed for reproducibility
par(bg="black") # set background to black

plot(x = clp, y = net.graph, edge.width = 0.1,  mark.border=NA,  vertex.label.cex=0.5, 
      vertex.label.color="white", edge.arrow.size=0.1, vertex.size = 6, edge.color = "Grey",
     col = communityColors, layout=layout_with_fr)

title(main = list("Animal's association network", cex = 1, col = "white", font = 2))

legend(x =-1.8 , y = 0.5, legend=levels(factor(names, levels=names)),
       col = unique(communityColors), bty = "n", pch=20 , pt.cex = 1.0, cex = 0.6, 
       text.col=unique(communityColors), horiz = FALSE, inset = c(0.1, 0.1), text.font =2)
dev.off()

# GRAPH 4
# Let's add leyends to the groups by picking the first animal on each cluster
# Plot with labels, clustered by community, and nodes scaled by weight (n# of links)

set.seed(28) # set seed for reproducibility
par(bg="black") # set background to black
plot(x = clp, y = net.graph, edge.width = 0.1,  mark.border=NA, vertex.label.cex=0.5, 
      vertex.label.color="white", edge.arrow.size=0.1, edge.color = "Grey",
     col = communityColors, layout=layout_with_fr)

title(main = list("Animal's association network", cex = 1, col = "white", font = 2))

legend(x =-1.6 , y = 0.5, legend=levels(factor(names, levels=names)),
       col = unique(communityColors), bty = "n", pch=20 , pt.cex = 1.5, cex = 0.7, 
       text.col=unique(communityColors), horiz = FALSE, inset = c(0.1, 0.1), text.font =2)

# Let's plot out GRAPH 4
tiff("plot_weighted_nodes__black_b_leg.tiff", width =6, height = 5, units = 'in', res = 300, compression = 'none')

set.seed(28) # set seed for reproducibility
par(bg="black") # set background to black

plot(x = clp, y = net.graph, edge.width = 0.1,  mark.border=NA, vertex.label.cex=0.5, 
      vertex.label.color="white", edge.arrow.size=0.1, edge.color = "Grey",
     col = communityColors, layout=layout_with_fr)

title(main = list("Animal's association network", cex = 1, col = "white", font = 2))

legend(x =-1.8 , y = 0.5, legend=levels(factor(names, levels=names)),
       col = unique(communityColors), bty = "n", pch=20 , pt.cex = 1.0, cex = 0.6, 
       text.col=unique(communityColors), horiz = FALSE, inset = c(0.1, 0.1), text.font =2)
dev.off()