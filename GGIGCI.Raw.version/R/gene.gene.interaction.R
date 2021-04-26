#!/usr/bin/env Rscript

#setwd("~/Desktop/data_from_HG/my_own_package/GGIGCI/")

gene.gene.interaction <- function(data, n.times = 1000, selected.gene.list=NULL, gene.level = c("gene.level", "default", "minP")){

  if (!inherits(data, "GGIGCI.data")) {
    stop("data must be one of the class 'GGIGCI.data'(the output of the function impute.res).")
  }
  snp.info <- data$SnpMatrix
  gene.info <- data$gene.info
  sample.info <- data$sample.info
  gene.level <- match.arg(gene.level)
  #phen <- sample.info[,6]
  if (is.null(selected.gene.list)) {
    selected.gene.list <- unique(gene.info[, 2])
  }
  else if (any(!selected.gene.list%in%gene.info[,2])) {
    stop("invaild gene selected.")
  }
  #else {
  gene.info <- gene.info[which(gene.info[, 2]%in%selected.gene.list), ]
  #}
  if (any(sample.info[,2]!=rownames(snp.info))){
    warning("miss match between sample info matrix and SnpMatrix.")
  }
  phen <- sample.info[,6]
  p.val.matrix <- diag(0, nrow = length(selected.gene.list))
  colnames(p.val.matrix) <- selected.gene.list
  rownames(p.val.matrix) <- selected.gene.list
  s.val.matrix <- diag(0, nrow = length(selected.gene.list))
  colnames(s.val.matrix) <- selected.gene.list
  rownames(s.val.matrix) <- selected.gene.list

  interaction.pairs <- combn(selected.gene.list, 2)
  gene.tabe <- data.frame(Genenames = selected.gene.list, Snps = table(gene.info[,2])[selected.gene.list])[,-2]
  ggi.ud <- data.frame(from = interaction.pairs[1, ], to = interaction.pairs[2, ],
                       P.val = rep(0, length(interaction.pairs[1, ])), S.val = rep(0, length(interaction.pairs[1, ])))
  nc.i <- ncol(interaction.pairs)

  bar <- txtProgressBar(0, nc.i, char = ">", style = 3)
  for(n in 1:nc.i){
    gene1 <- interaction.pairs[1, n]
    gene2 <- interaction.pairs[2, n]
    #print(gene1)
    #print(gene2)
    snp.gene1 <- gene.info[which(gene.info[, 2]==gene1), 3]
    snp.gene2 <- gene.info[which(gene.info[, 2]==gene2), 3]
    snp.info.gene1 <- snp.info[, colnames(snp.info)%in%snp.gene1]
    snp.info.gene2 <- snp.info[, colnames(snp.info)%in%snp.gene2]
    #temp.result <- ifelse(gene.level, GBIGM.gene.level(phen = phen, gene1 = snp.info.gene1, gene2 = snp.info.gene2, n.times = n.times),
                          #GBIGM(phen = phen, gene1 = snp.info.gene1, gene2 = snp.info.gene2, n.times = n.times))
    if (gene.level == "gene.level") {
      temp.result <- GBIGM.gene.level(phen = phen, gene1 = snp.info.gene1, gene2 = snp.info.gene2, n.times = n.times)
    }
    else if (gene.level == "minP") {
      #snp.info.gene1 <- fex(snp.info.gene1, gene.info[which(gene.info[, 2]==gene1), ])
      #snp.info.gene2 <- fex(snp.info.gene2, gene.info[which(gene.info[, 2]==gene2), ])
      division.gene1 <- cluster.gene(snp.info.gene1)
      division.gene2 <- cluster.gene(snp.info.gene2)

      pair.start <- expand.grid(division.gene1$start, division.gene2$start)
      pair.end <- expand.grid(division.gene1$end, division.gene2$end)

      temp.pair <- data.frame(start.gene1 = pair.start[,1], end.gene1 = pair.end[, 1], start.gene2 = pair.start[,2], end.gene2 = pair.end[, 2])
      p.list <- rep(NA, nrow(temp.pair))
      s.list <- rep(NA, nrow(temp.pair))

      for(i in 1:nrow(temp.pair)){
        sub.gene1 <- (temp.pair$start.gene1[i]):(temp.pair$end.gene1[i])
        sub.gene2 <- (temp.pair$start.gene2[i]):(temp.pair$end.gene2[i])
        sub.p.s <- GBIGM(phen = phen, snp.info.gene1[, sub.gene1], snp.info.gene2[, sub.gene2], n.times = n.times)
        p.list[i] <- sub.p.s[1]
        s.list[i] <- sub.p.s[2]
      }
      
      p.list <- p.adjust(p.list, "BH")
      temp.result <- c(min(p.list), s.list[which.min(p.list)])
    }
    else if (gene.level == "default") {
      temp.result <- GBIGM(phen = phen, gene1 = snp.info.gene1, gene2 = snp.info.gene2, n.times = n.times)
      #temp.result <- GBIGM.test(Y = phen, G1 = snp.info.gene1, G2 = snp.info.gene2, n.perm = n.times)
    }
    p.val.matrix[gene1, gene2] <- p.val.matrix[gene2, gene1] <- temp.result[1]
    #p.val.matrix[gene1, gene2] <- p.val.matrix[gene2, gene1] <- temp.result$p.val
    s.val.matrix[gene1, gene2] <- s.val.matrix[gene2, gene1] <- temp.result[2]
    #s.val.matrix[gene1, gene2] <- s.val.matrix[gene2, gene1] <- temp.result$statistic
    #ggi.ud[n, 3:4] <- c(temp.result$p.val, temp.result$statistic)
    ggi.ud[n, 3:4] <- temp.result

    setTxtProgressBar(bar, n)
  }

  message("\n")
  snp.matrix <- snp.info[, colnames(snp.info)%in%gene.info[,3]]
  out <- list(SnpMatrix = snp.matrix, gene.info = gene.info, sample.info = sample.info,
              GGI.ud = ggi.ud, GENE.table = gene.tabe, P.matrix = p.val.matrix, S.matrix = s.val.matrix)
  class(out) <- c("GGIUD", "list")
  return(out)
}

cluster.gene <- function(gene){
  if (ncol(gene)>20){
    distance.gene <- snpStats::ld(gene, gene, stats = "R.squared")
    distance.gene <- as.dist(1 - distance.gene)
    clust.tree.gene <- rioja::chclust(distance.gene)
    k.gene <- cutree(clust.tree.gene, k = 1:(ncol(gene)-19))
    max.gene <- sapply(1:(ncol(gene)-19),FUN=function(i){return(max(table(as.factor(k.gene[,i]))))})
    id.gene <- which(max.gene<=20)[1]
  }
  else if (ncol(gene)==1){
    return(list(start = c(1), end = c(1)))
  }
  else {
    k.gene <- matrix(rep(1, ncol(gene)), ncol = 1)
    rownames(k.gene) <- colnames(gene)
    max.gene <- ncol(gene)
    id.gene <- 1
  }
  #print(sapply(1:(length(k.gene[,id.gene])-1),FUN=function(i){return(k.gene[,id.gene][i]!=k.gene[,id.gene][i+1])}))

  division.gene.start <- c(1,1+as.numeric(which(sapply(1:(length(k.gene[,id.gene])-1),FUN=function(i){return(k.gene[,id.gene][i]!=k.gene[,id.gene][i+1])}))))
  division.gene.end <- c(as.numeric(which(sapply(1:(length(k.gene[,id.gene])-1),FUN=function(i){return(k.gene[,id.gene][i]!=k.gene[,id.gene][i+1])}))),ncol(gene))

  return(list(start = division.gene.start, end = division.gene.end))
}
