#load("~/Desktop/history.rdata/first.data.res.after.impute20201123.rda")
#gguid.gene.level.multi.loci <- gene.gene.interaction(temp.im, gene.level = T, selected.gene.list = gene.multi.loci)
#gguid.normal.multi.loci <- gene.gene.interaction(temp.im, gene.level = F, selected.gene.list = gene.multi.loci)
#save.image(file = "~/Desktop/history.rdata/my.image.20201204.RData")

#gguid.ra.gbigm <- gene.gene.interaction(res.ra.after.impute)
#evex.official.symbol <- data.frame(entrezgene_id = evex.id.list, official_symbol = rep(NA, length(evex.id.list)))
#for(id in evex.id.list) {
#  if (dim(evex.symbol[evex.symbol[, 1]==id & evex.symbol[,2]=="official_symbol", ])[1]==0) {
#    evex.official.symbol[evex.official.symbol[, 1]==id, 2] <- paste0("EVEX_", as.character(id))
#  }
#  else {
#    temp <- evex.symbol[evex.symbol[, 1]==id & evex.symbol[,2]=="official_symbol", ]
#    evex.official.symbol[evex.official.symbol[, 1]==id, 2] <- temp[1, 3]
#  }
#}

#gguid.all.add.minP.adjust <- gene.gene.interaction(res, gene.level = "minP")

#load("~/Desktop/history.rdata/protein.network.with.3Dcoordinate.RData")
#
#gene.list <- read.table("~/Desktop/data_from_HG/string/first_data_17_gene/genename17.txt", 
#                        header = F)
#require(igraph)
#string.graph <- graph_from_data_frame(d = protein.links, directed = F, 
#                                      vertices = all.protein)
#gene.combn <- as.data.frame(t(combn(gene.list[,1], 2)))
#colnames(gene.combn) <- c("protein1", "protein2")
#gene.combn$step <- NA
#gene.combn$distance <- NA
#
#for (i in 1:dim(gene.combn)[1]) {
#  gene1 <- gene.combn[i, 1]
#  gene2 <- gene.combn[i, 2]
#  message(paste0("the shortest routes between ", gene1, " and ", gene2))
#  protein1 <- all.protein$ensp.id[which(all.protein$gene.symbol==gene1)]
#  protein2 <- all.protein$ensp.id[which(all.protein$gene.symbol==gene2)]
#  #print(paste0(protein1, " and ", protein2))
#  #path1to2 <- shortest_paths(string.graph, from = protein1, to = protein2, 
#  #mode = "all", weights = E(string.graph)$distance)$vpath
#  #print(path1to2)
#  #for (path in path1to2){
#  #  idx <- as.vector(path)
#  #  route <- names(V(string.graph)[idx])
#  #  #print(route)
#  #  route <- all.protein$gene.symbol[which(all.protein$ensp.id%in%route)]
#  #  print(route)
#  #}
#  
#  path1to2 <- shortest_paths(string.graph, from = protein1, to = protein2, 
#                             mode = "all", weights = NA, output = "both")$epath[[1]]
#  idx <- as.vector(path1to2)
#  gene.combn$step[i] <- length(idx)
#  gene.combn$distance[i] <- sum(protein.links$distance[idx])
# 
#  print(protein.links[idx, ])
#  message("distance: ", sum(protein.links$distance[idx]))
#} 

#print(gene.combn)

#gguid.normal.default <- gene.gene.interaction(temp.im, gene.level = "default")
#gguid.third.c <- gene.gene.interaction(res.third, gene.level = "default")
#gguid.third.adminp <- gene.gene.interaction(res.third, gene.level = "minP")
pair.third.496 <- as.data.frame(t(combn(gene.third.in.string$preferred_name, 2)))
pair.third.496$distance.1 <- NA
pair.third.496$distance.2 <- NA
pair.third.496$distance.3 <- NA
colnames(pair.third.496) <- c("from", "to", "distance.step", "distance.raw", "distance.distance")
require(igraph)
for(i in 1:nrow(pair.third.496)) {
  message(i)
  gene1 <- pair.third.496[i, 1]
  gene2 <- pair.third.496[i, 2]
  gene1 <- gene.third.in.string[gene.third.in.string$preferred_name==gene1, 1]
  gene2 <- gene.third.in.string[gene.third.in.string$preferred_name==gene2, 1]
  e.path.shortest <- shortest_paths(network.string, from = gene1, to = gene2, mode = "all", weights = NA, output = "epath")$epath[[1]]
  distance <- sum(e.path.shortest$distance)
  pair.third.496$distance.step[i] <- distance
  e.path.shortest <- shortest_paths(network.string, from = gene1, to = gene2, mode = "all", output = "epath")$epath[[1]]
  distance <- sum(e.path.shortest$distance)
  pair.third.496$distance.raw[i] <- distance
  e.path.shortest <- shortest_paths(network.string, from = gene1, to = gene2, mode = "all", weights = protein.links$distance, output = "epath")$epath[[1]]
  distance <- sum(e.path.shortest$distance)
  pair.third.496$distance.distance[i] <- distance
}

save(pair.third.496, file = "~/Desktop/data_from_HG/data_third_20210123/temp.pair.rda")

for(gene in gene.multi.loci){
  message(gene)
  snp.gene <- gene.info[gene.info$Genenames==gene, 3]
  snp.gene <- temp.im$SnpMatrix[, snp.gene]
  snp.gene <- as(snp.gene, "numeric")
  message(ncol(snp.gene))
  message(ncol(fex(snp.gene)))
}
