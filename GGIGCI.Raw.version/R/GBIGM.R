#!/usr/bin/env Rscript

#setwd("~/Desktop/data_from_HG/my_own_package/GGIGCI/")

GBIGM <- function(phen, gene1, gene2, n.times = 1000) {
  x1 <- as(gene1, "numeric")
  x2 <- as(gene2, "numeric")
  #selfhead(x1)
  #selfhead(x2)
  if (nrow(x1)!=nrow(x2)){
    stop("The two SnpMatrix must have the same row numbers. ")
  }
  else if (length(phen)!=nrow(x1)){
    stop("The length of phen must be equal to the SnpMatrix's row.number.")
  }
  num.sample <- length(phen)

  g1 <- apply(x1, 1, paste, collapse = " ")
  g2 <- apply(x2, 1, paste, collapse = " ")
  g1.2 <- apply(cbind(x1, x2), 1, paste, collapse = " ")

  gx1 <- apply(cbind(phen, x1), 1, paste, collapse = " ")
  gx2 <- apply(cbind(phen, x2), 1, paste, collapse = " ")
  gx1.2 <- apply(cbind(phen, cbind(x1, x2)), 1, paste, collapse = " ")
  #print(table(gx1))
  #print(table(gx2))
  #print(table(gx1.2))
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
    gs1 <- apply(cbind(D[[i]], x1), 1, paste, collapse = " ")
    #gs2 <- gx2[D[[i]]]
    gs2 <- apply(cbind(D[[i]], x2), 1, paste, collapse = " ")
    #gs1.2 <- gx1.2[D[[i]]]
    gs1.2 <- apply(cbind(D[[i]], cbind(x1, x2)), 1, paste, collapse = " ")

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
