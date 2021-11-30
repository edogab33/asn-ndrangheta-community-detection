#!/usr/bin/Rscript
## Author: Edoardo Gabrielli
## Date: 30 Nov 2021

library(igraph)
library(networkdata)
library(RColorBrewer)
library(vtable)

shrinkGraph <- function(g, membership) {
  ## Author: Alexander Stivala
  V(g)$membership <- membership
  gshrink <- contract.vertices(g, V(g)$membership,
                               vertex.attr.comb=list(membership='first',
                                                     lat='mean',
                                                     longi='mean',
                                                     discharges='sum',
                                                     states='first',
                                                     'ignore'))
  gshrink <- simplify(gshrink, remove.multiple=TRUE, remove.loops=TRUE,
                      edge.attr.comb='sum')
  return(gshrink)
}

assign_colors <- function(community) {
  ## Assign colors to nodes based on their community
  ## Author: Jelena Čuklina
  ## Source: https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
  n <- max(membership(community))
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  colors <- sample(col_vector, n)
  return(colors)
}


## Load 'Ndrangheta dataset
g <- covert_33
shape = ifelse(V(covert_33)$type, "square", "circle")

plot(g, vertex.shape=shape, vertex.color=V(g)$type, vertex.label = NA, vertex.size=3)

## Basic statistics
vcount(g)
summary(degree(g, V(g)$type==FALSE))
summary(degree(g, V(g)$type==TRUE))
graph.density(g)
hist(degree(g, V(g)$type==FALSE))
hist(degree(g, V(g)$type==TRUE))

## Projection to one-mode network
onemode = bipartite.projection(g, multiplicity = TRUE, which=FALSE)
plot(onemode, vertex.shape=shape, vertex.color=V(g)$type, vertex.label = NA, vertex.size=3)

## Basic statistics of one-mode network
vcount(onemode)
graph.density(onemode)
hist(degree(onemode))
summary(degree(onemode))

## Calculate clusters based on Walktrap and Louvain
cl <- cluster_walktrap(onemode, steps = 6)
print(modularity(cl))

cl2 <- cluster_louvain(onemode)
print(modularity(cl2))

# Plot walktrap clusters and Louvain ones
colors_cl <- assign_colors(cl)
l <- layout_nicely(onemode)
plot(onemode, vertex.shape=shape, vertex.color=colors_cl[membership(cl)], vertex.label = NA, vertex.size=3, layout = l)

colors_cl2 <- assign_colors(cl2)
plot(onemode, vertex.shape=shape, vertex.color=colors_cl2[membership(cl2)], vertex.label = NA, vertex.size=3, layout = l)

## Compare the two communities based on NMI
compare(cl, cl2, method = "nmi")

## Number of communities and community sizes
length(cl)
length(cl2)

sizes(cl)
sizes(cl2)

## Distribution of the number of nodes in each community
h1 <- hist(sizes(cl))
h2 <- hist(sizes(cl2))
plot(h1, col=rgb(0,0,1,1/4), xlim=c(0,60))
plot(h2, col=rgb(1,0,0,1/4), xlim=c(0,60), add=T)


### More experiments
## Shrink
gs <- shrinkGraph(onemode, membership(cl))
plot(gs)

## Take the subgraph of the biggest community
sub = induced.subgraph(onemode, is.element(V(onemode)$name, cl2[[which.max(sizes(cl2))]]))
plot(sub, vertex.shape=shape, vertex.label = NA, vertex.size=3)

## Cluster the subgraph
cl3 <- cluster_louvain(sub)
print(modularity(cl3))
colors_cl <- assign_colors(cl3)
lay <- layout_nicely(sub)
plot(sub, vertex.shape=shape, vertex.color=colors_cl[membership(cl)], vertex.label = NA, vertex.size=3, layout = lay)
