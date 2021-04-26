#!/usr/bin/env Rscript

#setwd("~/Desktop/data_from_HG/my_own_package/GGIGCIraw/")

plot.GGC.data <- function(ggc, color.method = c("snp", "out"), color.threshold = 5, shape.method = c("snp", "out"), shape.threshold = 5){
  #require(igraph)
  color.method <- match.arg(color.method)
  shape.method <- match.arg(shape.method)

  graph.data <- ggc$gg.causal
  graph.vertices <- ggc$gene.table
  ngenes <- dim(graph.vertices)[1]

  graph.vertices$out.edge <- rep(0, ngenes)
  for (gene in graph.vertices[,1]){
    graph.vertices[which(graph.vertices[,1]==gene),]$out.edge <- sum(graph.data[,1]==gene)
  }

  graph.vertices$color <- rep("", ngenes)
  if (color.method=="snp"){
    #graph.vertices[which(graph.vertices$Snps.Freq>color.threshold),]$color <- "red"
    #graph.vertices[which(graph.vertices$Snps.Freq<=color.threshold),]$color <- "blue"
    graph.vertices$color <- ifelse(graph.vertices$Snps.Freq>color.threshold, "red", "yellow")
  }
  else if (color.method=="out"){
    #graph.vertices[which(graph.vertices$out.edge>color.threshold),]$color <- "red"
    #graph.vertices[which(graph.vertices$out.edge<=color.threshold),]$color <- "blue"
    graph.vertices$color <- ifelse(graph.vertices$out.edge>color.threshold, "red", "yellow")
  }

  graph.vertices$shape <- rep("", ngenes)
  if (shape.method=="snp"){
    #graph.vertices[which(graph.vertices$Snps.Freq>shape.threshold), ]$shape <- "rectangle"
    #graph.vertices[which(graph.vertices$Snps.Freq<=shape.threshold),]$shape <- "pie"
    graph.vertices$shape <- ifelse(graph.vertices$Snps.Freq>shape.threshold, "rectangle", "circle")
  }
  else if (shape.method=="out"){
    #graph.vertices[which(graph.vertices$out.edge>shape.threshold), ]$shape <- "rectangle"
    #graph.vertices[which(graph.vertices$out.edge<=shape.threshold),]$shape <- "pie"
    graph.vertices$shape <- ifelse(graph.vertices$out.edge>shape.threshold, "rectangle", "circle")
  }

  print(graph.vertices)
  print(graph.data)

  network.ggc <- igraph::graph_from_data_frame(graph.data, directed = T, vertices = graph.vertices)
  plot(network.ggc, vertex.shapes = V(network.ggc)$shape, vertex.color = V(network.ggc)$color)
}
