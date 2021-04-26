#!/usr/bin/env Rscript

#setwd("~/Desktop/data_from_HG/my_own_package/GGIGCI/")

gene.gene.causal <- function(data.ggiud, n.times = 1000, sample.size = 100,
                             ref = c("Gaussian", "uniform"), est = c("entropy", "intergral"),
                             interaction.p = 0.05, fex = T, pos.use = T) {
  #if (!inherits(data.ggigci, "GGIGCI.data")){
  #  stop("data must be one of the class 'GGIGCI.data'(the output of the function impute.res).")
  #}
  if (!inherits(data.ggiud, "GGIUD")){
    stop("data must be one of the class 'GGIUD'(the output of the function gene.gene.interaction).")
  }

  ref <- match.arg(ref)
  est <- match.arg(est)
  #print(paste0(ref, "xiix"))
  #print(est)
  gene.table <- data.ggiud$GENE.table
  snp.info <- data.ggiud$SnpMatrix
  gene.info <- data.ggiud$gene.info
  sample.info <- data.ggiud$sample.info
  phen <- as(sample.info[,6], "numeric")
  selected.gene.pairs <- data.ggiud$GGI.ud[which(data.ggiud$GGI.ud$P.val<=interaction.p), ]
  gene.f.table <- data.frame(gene1 = selected.gene.pairs[, 1], gene2 = selected.gene.pairs[, 2],
                             f1 = rep(NA, length(selected.gene.pairs[,1])), f2 = rep(NA, length(selected.gene.pairs[, 2])),
                             P.value = rep(0, length(selected.gene.pairs[, 1])), causal.direction = rep(0, length(selected.gene.pairs[, 1])))
  gene.causal.table <- data.frame(from = selected.gene.pairs[, 1], to = selected.gene.pairs[, 2],
                                  F.value = rep(0, length(selected.gene.pairs[, 1])))
  n.pairs <- dim(selected.gene.pairs)[1]
  for(i in seq_len(n.pairs)){
    #print(selected.gene.pairs[i,])
    gene1 <- selected.gene.pairs[i, 1]
    gene2 <- selected.gene.pairs[i, 2]
    snp.gene1 <- gene.info[which(gene.info[,2]==gene1), 3]
    snp.gene2 <- gene.info[which(gene.info[,2]==gene2), 3]
    snp.info.gene1 <- as(snp.info[, colnames(snp.info)%in%snp.gene1], "numeric")
    snp.info.gene2 <- as(snp.info[, colnames(snp.info)%in%snp.gene2], "numeric")

    if (fex){
      #x.gene1 <- fex(cbind(phen, snp.info.gene1))
      #x.gene2 <- fex(cbind(phen, snp.info.gene2))
      if(pos.use){
        x.gene1 <- fex(snp.info.gene1, gene.info[which(gene.info[,2]==gene1), ])
        x.gene2 <- fex(snp.info.gene2, gene.info[which(gene.info[,2]==gene2), ])
      }
      else {
        x.gene1 <- fex(snp.info.gene1)
        x.gene2 <- fex(snp.info.gene2)
      }

      F.matrix <- diag(0, nrow = dim(x.gene1)[2], ncol = dim(x.gene2)[2])
      rownames(F.matrix) <- colnames(x.gene1)
      colnames(F.matrix) <- colnames(x.gene2)
      for(n in seq_len(dim(F.matrix)[1])) {
        for(m in seq_len(dim(F.matrix)[2])) {
          pairs.data <- data.frame(gene1 = x.gene1[, n], gene2 = x.gene2[, m])
          #fex.temp.result <- IGCIfor2value(pairs.data, times = n.times, sample.size = sample.size, refMeasure = ref, estimator = est)
          F.matrix[n, m] <- igci(pairs.data[,1], pairs.data[, 2], refMeasure = ref, estimator = est)
        }
      }
      F1 <- sum(F.matrix<0)
      F2 <- sum(F.matrix>0)
      gene.f.table[i, 3:5] <- c(F1, F2, NA)
      if (F1>F2){
        gene.f.table[i, 5] <- F2/(F1+F2)
        gene.f.table[i, 6] <- "gene1->gene2"
      }
      else if (F1<F2) {
        gene.f.table[i, 5] <- F1/(F1+F2)
        gene.f.table[i, 6] <- "gene2->gene1"
      }
    }
    else {
      x.gene1 <- apply(cbind(phen, 2*base3to10(snp.info.gene1)), 1, sum)
      x.gene2 <- apply(cbind(phen, 2*base3to10(snp.info.gene2)), 1, sum)
      pairs.data <- data.frame(gene1 = x.gene1, gene2 = x.gene2)
      temp.result <- IGCIfor2value(pairs.data, times = n.times, sample.size = sample.size, refMeasure = ref, estimator = est)

      gene.f.table[i, 3:6] <- c(temp.result$f.mean[1], temp.result$f.mean[2], temp.result$p.value, temp.result$causal.direction)
    }
    if (gene.f.table[i, 6]=="gene1->gene2"){
      gene.causal.table[i, ] <- c(gene.f.table[i, 1], gene.f.table[i, 2], gene.f.table[i, 3])
    }
    else if (gene.f.table[i, 6]=="gene2->gene1") {
      gene.causal.table[i, ] <- c(gene.f.table[i, 2], gene.f.table[i, 1], gene.f.table[i, 4])
    }
  }
  out <- list(SnpMatrix = snp.info, gene.info = gene.info, sample.info = sample.info,
              gg.causal = gene.causal.table, gene.table = gene.table, gg.f = gene.f.table)
  if (fex){
    out$method <- "Fexpand"
  }
  else {
    out$method <- "base3to10"
  }

  class(out) <- c("GGC.data", "list")
  return(out)
}

#fex <- function(gene, pos){
#
#  return(apply(gene, 1, sum))
#}
fex <- function(gene, pos=NULL) {
  if (is.null(pos)) {
    pos.v <- (0:(dim(gene)[2]-1))/(dim(gene)[2]-1)
    names(pos.v) <- colnames(gene)
    gene.name <- "gene"
  }
  else {
    pos.v <- as(pos$Position, "numeric")
    names(pos.v) <- pos$SNPnames
    gene.name <- pos$Genenames[1]
  }
  nsnps <- length(pos.v)
  if (nsnps>1){
    idx <- order(pos.v)
    gene <- gene[, idx]
    pos.v <- pos.v[idx]
    pos.v <- (pos.v - pos.v[1])/(pos.v[nsnps]- pos.v[1])
    #print(pos.v)

    eigenval <- prcomp(gene)$sd^2
    eigen.sum <- sum(eigenval)
    tmp <- 0
    n.of.basis <- 0

    for (i in 1:length(eigenval)){
      tmp <- eigenval[i]+tmp
      n.of.basis <- i
      if (tmp >= 0.8*eigen.sum){
        break
      }
    }
    n.of.basis <- floor(n.of.basis/2)*2+1
    #print(n.of.basis)
    frange <- c(pos.v[1], pos.v[length(pos.v)])
    fbasis <- fda::create.fourier.basis(frange, nbasis = n.of.basis)
    phi <- fda::eval.basis(pos.v, fbasis)
    res <- t(MASS::ginv(t(phi)%*%phi)%*%t(phi)%*%t(gene))
  }
  else {
    res <- gene
  }
  col.num <- dim(res)[2]
  #gene.name <- pos$Genenames[1]
  #print(gene.name)
  colnames(res) <- paste(gene.name, seq_len(col.num), sep = ".")

  return(res)
}

base3to10 <- function(gene){
  n.SNP<-ncol(gene)
  long<-nrow(gene)
  tab<-matrix(NA,ncol=n.SNP,nrow=long)
  for (i in seq_len(n.SNP)){
    tab[,i]<-3^(i-1)*gene[,i]
  }
  b10<-apply(tab,1,sum)
  return(b10)
}
