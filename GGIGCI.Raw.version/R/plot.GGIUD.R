#!/usr/bin/env Rscript

#setwd("~/Desktop/data_from_HG/my_own_package/GGIGCIraw/")

plot.GGIUD <- function(ggiud, p.threshold = 0.05, color.method = c("snp", "degree"), color.threshold = 5, shape.method = c("snp", "degree"), shape.threshold = 5) {
  require(igraph)

  color.method <- match.arg(color.method)
  shape.method <- match.arg(shape.method)

  graph.data <- ggiud$GGI.ud[ggiud$GGI.ud$P.val<p.threshold, ]
  graph.vertices <- ggiud$GENE.table
  ngenes <- dim(graph.vertices)[1]

  graph.vertices$degree <- rep(0, ngenes)
  for (gene in graph.vertices[,1]){
    graph.vertices[which(graph.vertices[,1]==gene),]$degree <- sum(graph.data[,1]==gene) + sum(graph.data[,2]==gene)
  }
  graph.vertices$color <- rep("", ngenes)
  if (color.method=="snp"){
    graph.vertices$color <- ifelse(graph.vertices$Snps.Freq > color.threshold, "red", "yellow")
  }
  else if (color.method=="degree"){
    graph.vertices$color <- ifelse(graph.vertices$degree > color.threshold, "red", "yellow")
  }

  graph.vertices$shape <- rep("", ngenes)
  if (shape.method=="snp"){
    graph.vertices$shape <- ifelse(graph.vertices$Snps.Freq>shape.threshold, "rectangle", "circle")
  }
  else if (shape.method=="degree"){
    graph.vertices$shape <- ifelse(graph.vertices$degree > shape.threshold, "rectangle", "circle")
  }

  print(graph.vertices)
  print(graph.data)

  network.ggiud <- igraph::graph_from_data_frame(graph.data, directed = F, vertices = graph.vertices)
  plot(network.ggiud, vertex.shapes = V(network.ggiud)$shape, vertex.color = V(network.ggiud)$color)
}
