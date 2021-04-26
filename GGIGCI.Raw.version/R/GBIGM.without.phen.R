#!/usr/bin/env Rscript

#setwd("~/Desktop/data_from_HG/my_own_package/GGIGCI/")

GBIGM.without.phen <- function(gene1, gene2, n.times = 1000){
  x1 <- as(gene1, "numeric")
  x2 <- as(gene2, "numeric")
  #selfhead(x1)
  #selfhead(x2)
  if (nrow(x1)!=nrow(x2)){
    stop("The two SnpMatrix must have the same row numbers. ")
  }
  num.sample <- nrow(gene1)

  gx1 <- apply(x1, 1, paste, collapse = " ")
  gx2 <- apply(x2, 1, paste, collapse = " ")
  gx1.2 <- apply(cbind(x1, x2), 1, paste, collapse = " ")
  #print(table(gx1))
  #print(table(gx2))
  HG1 <- entropy.vec(gx1)
  HG2 <- entropy.vec(gx2)
  HG1.2 <- entropy.vec(gx1.2)

  if (min(HG1, HG2)==0){
    Delta1.2init <- min(HG1, HG2)-HG1.2
  }
  else {
    Delta1.2init <- (min(HG1, HG2)-HG1.2)/min(HG1, HG2)
  }

  D <- list()
  Delta1.2 <- rep(NA, n.times)
  for (i in seq_len(n.times)){
    D[[i]] <- sample(1:num.sample, num.sample)
    gs1 <- gx1[D[[i]]]
    gs2 <- gx2[D[[i]]]
    gs1.2 <- gx1.2[D[[i]]]
    HGs1 <- entropy.vec(gs1)
    HGs2 <- entropy.vec(gs2)
    HGs1.2 <- entropy.vec(gs1.2)
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

  return(list(num.sample, c(pval, stat), Delta1.2))
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
