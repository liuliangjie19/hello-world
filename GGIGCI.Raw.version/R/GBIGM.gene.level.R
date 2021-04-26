#!/usr/bin/env Rscript

#setwd()
###when the number of snp loci on the your gene is too small, the gene.level.method can be considered.
GBIGM.gene.level <- function(phen, gene1, gene2, n.times = 1000) {
  
  x1 <- as(gene1, "numeric")
  x2 <- as(gene2, "numeric")
  xt1 <- apply(x1, 2, correct.genotype)
  xt2 <- apply(x2, 2, correct.genotype)
  
  if (nrow(x1)!=nrow(x2)){
    stop("Error 1.")
  }
  else if (length(phen)!=nrow(x1)){
    stop("Error 2")
  }
  num.sample <- length(phen)
  
  g1 <- apply(xt1, 1, is.mutation)
  g2 <- apply(xt2, 1, is.mutation)
  #g1.2 <- apply(cbind(xt1, xt2), 1, is.mutation)
  g1.2 <- apply(cbind(g1, g2), 1, paste, collapse = " ")
  
  gx1 <- apply(cbind(phen, g1), 1, paste, collapse = " ")
  gx2 <- apply(cbind(phen, g2), 1, paste, collapse = " ")
  gx1.2 <- apply(cbind(phen, g1.2), 1, paste, collapse = " ")
  
  HG1 <- entropy.vec(g1)
  HG2 <- entropy.vec(g2)
  HG1.2 <- entropy.vec(g1.2)
  
  HGx1 <- entropy.vec(gx1)-HG1
  HGx2 <- entropy.vec(gx2)-HG2
  HGx1.2 <- entropy.vec(gx1.2)-HG1.2
  
  if (min(HGx1, HGx2)==0){
    Delta1.2init <- min(HGx1, HGx2)-HGx1.2
  }
  else {
    Delta1.2init <- (min(HGx1, HGx2)-HGx1.2)/min(HGx1, HGx2)
  }
  
  D <- list()
  for(i in seq_len(n.times)){
    D[[i]] <- phen[sample(1:num.sample, num.sample)]
  }
  Delta1.2 <- rep(NA, n.times)
  for (i in seq_len(n.times)){
    #D[[i]] <- sample(1:num.sample, num.sample)
    #gs1 <- gx1[D[[i]]]
    gs1 <- apply(cbind(D[[i]], g1), 1, paste, collapse = " ")
    #gs2 <- gx2[D[[i]]]
    gs2 <- apply(cbind(D[[1]], g2), 1, paste, collapse = " ")
    #gs1.2 <- gx1.2[D[[i]]]
    gs1.2 <- apply(cbind(D[[i]], g1.2), 1, paste, collapse = " ")
    
    HGs1 <- entropy.vec(gs1)-HG1
    HGs2 <- entropy.vec(gs2)-HG2
    HGs1.2 <- entropy.vec(gs1.2)-HG1.2
    if (min(HGs1, HGs2)==0){
      Delta1.2tmp <- min(HGs1, HGs2)-HGs1.2
    }
    else {
      Delta1.2tmp <- (min(HGs1, HGs2)-HGs1.2)/min(HGs1, HGs2)
    }
    Delta1.2[i] <- Delta1.2tmp
  }
  pval <- sum(Delta1.2>=Delta1.2init)/n.times
  
  stat <- Delta1.2init
  
  return(c(pval, stat))
}

correct.genotype <- function(snp.x){
  w <- "Wildtype"
  m <- "Mutation"
  
  #snp.table <- table(snp.x)
  l0 <- length(snp.x[which(snp.x==0)])
  l2 <- length(snp.x[which(snp.x==2)])
  #if (snp.table[[1]]<snp.table[[2]]){
  if (l0<l2){
    snp.t <- ifelse(snp.x==2, w, m)
  }
  else {
    snp.t <- ifelse(snp.x==0, w, m)
  }
  
  return(snp.t)
}

is.mutation <- function(xt){
  return(ifelse("Wildtype"%in%xt, 1, 0))
}

entropy.vec <- function(X){
  pourc<-percentofX(X)
  H<-(-sum(pourc*log(pourc)))
  return(H)
}

percentofX <- function(X){
  n.ind<-length(X)
  b10<-as.factor(X)
  vect<-tabulate(b10)
  pourc<-vect/n.ind
  return(pourc)
}
