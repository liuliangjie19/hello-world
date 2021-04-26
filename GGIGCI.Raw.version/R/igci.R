# R version of Information Geometry Causal Inference
#                _                           
# platform       x86_64-apple-darwin17.0     
# arch           x86_64                      
# os             darwin17.0                  
# system         x86_64, darwin17.0          
# status                                     
# major          4                           
# minor          0.0                         
# year           2020                        
# month          04                          
# day            24                          
# svn rev        78286                       
# language       R                           
# version.string R version 4.0.0 (2020-04-24)
# nickname       Arbor Day 

# P. Daniusis, D. Janzing, J. Mooij, J. Zscheischler, B. Steudel,
# K. Zhang, B. Scholkopf:  Inferring deterministic causal relations.
# Proceedings of the 26th Annual Conference on Uncertainty in Artificial 
# Intelligence (UAI-2010).  
# http://event.cwi.nl/uai2010/papers/UAI2010_0121.pdf
#######################################################################################
# Inputs:                                                                             #
#   x                       L x 1 observations of x                                   #
#   y                       L x 1 observations of y                                   #
# refMeasure      reference measure to use:                                           #
#   1: uniform                                                                        #
#   2: Gaussian                                                                       #
# estimator       estimator to use:                                                   #
#   1: entropy (eq. (12) in [1]),                                                     #
#   2: integral approximation (eq. (13) in [1]).                                      #
#                                                                                     #
# Outputs:                                                                            #
#   f < 0:          the method prefers the causal direction x -> y                    #
#   f > 0:          the method prefers the causal direction y -> x                    #
#######################################################################################


psi <- function(x) {
  dx <- D(expression(log(gamma(x))), "x")
  return(eval(dx))
}

igci <- function(x, y, refMeasure = "Gaussian", estimator = "entropy") {
  if (refMeasure=="Gaussian" || refMeasure=="g") {
    refMeasure <- 2
  } else if (refMeasure == "uniform" || refMeasure=="u") {
    refMeasure <- 1
  } else {
    stop("refMeasure must be Gaussian or uniform.\n")
  }
  if (estimator=="entropy" || estimator=="e") {
    estimator <- 1
  } else if (estimator == "intergral" || estimator=="i") {
    estimator <- 2
  } else {
    stop("estimator must be extropy or intergral.\n")
  }
  if (!(is.vector(x) & is.vector(y))) {
    stop("x and y must be vector!\n")
  }
  nx <- length(x)
  ny <- length(y)
  if (nx != ny) {
    stop("the length of x and y must be equal!\n")
  } 
  if (nx < 20) {
    stop("Not enough observations in this case!\n")
  }
  
  if (refMeasure == 1) {
    xi <- (x-min(x))/(max(x)-min(x))
    yi <- (y-min(y))/(max(y)-min(y))
  } else {
    xi <- (x-mean(x))/(sd(x))
    yi <- (y-mean(y))/(sd(y))
  }
  x.y <- data.frame(data.x = xi, data.y = yi)
  f <- 0
  if (estimator == 1) {
    x.sort <- sort(x.y$data.x)
    y.sort <- sort(x.y$data.y)
    hx <- 0
    for(i in 1:(nx-1)) {
      delta <- x.sort[i+1] - x.sort[i]
      if (delta) {
        hx <- hx + log(abs(delta))
      }
    }
    #hx <- hx/(nx-1)+psi(nx)-psi(1)
    hx <- hx/(nx-1) + digamma(nx) - digamma(1)
    
    hy <- 0
    for(i in 1:(ny-1)) {
      delta <- y.sort[i+1]- y.sort[i]
      if (delta){
        hy <- hy + log(abs(delta))
      }
    }
    #hy <- hy/(ny-1)+psi(ny)-psi(1)
    hy <- hy/(ny-1) + digamma(ny) - digamma(1)
    
    f <- hy - hx
    return(f)
    
  } else {
    a <- 0
    b <- 0
    
    indx <- sort(x.y$data.x, index.return = T)$ix
    indy <- sort(x.y$data.y, index.return = T)$ix
    x.sortx <- x.y[indx, ]$data.x
    y.sortx <- x.y[indx, ]$data.y
    x.sorty <- x.y[indy, ]$data.x
    y.sorty <- x.y[indy, ]$data.y
    for(i in 1:(nx-1)) {
      if (x.sortx[i+1]!=x.sortx[i] & y.sortx[i+1]!=y.sortx[i]){
        a <- a + log(abs(y.sortx[i+1]-y.sortx[i])/(x.sortx[i+1]-x.sortx[i]))
      }
      if (y.sorty[i+1]!=y.sorty[i] & x.sorty[i+1]!=x.sorty[i]){
        b <- b + log(abs(x.sorty[i+1]-x.sorty[i])/(y.sorty[i+1]-y.sorty[i]))
      }
    }
    
    f <- (a-b)/nx
    return(f)
  }
  
  #return(f)
}

IGCIfor2value <- function(data, sample.size = 100, times = 100, refMeasure = "g", estimator = "e") {
  if (dim(data)[2]!=2) {
    stop("More than 2 values!\n")
  }
  if (dim(data)[1]<100) {
    stop("less than 100 obervations!\n")
  }
  F1 <- 0
  F2 <- 0
  f1 <- 0
  f2 <- 0
  n <- dim(data)[1]
  for (i in 1:times) {
    ind.sample <- sample(n, sample.size, replace = T)
    data.sample <- data[ind.sample, ]
    x <- data.sample[,1]
    y <- data.sample[,2]
    f <- igci(x, y, refMeasure = refMeasure, estimator = estimator)
    if (f<0) {
      F1 <- F1 + 1
      f1 <- f1 + f
    } else if (f>0){
      F2 <- F2 + 1
      f2 <- f2 + f
    }else {
      next
    }
  }
  f.mean <- c(f1/F1, f2/F2)
  if (F1==0){
    f.mean[1] <- NA
  }
  if (F2==0){
    f.mean[2] <- NA
  }
  if (F1 > F2) {
    p.value <- 1 - F1/times
    causal <- paste0(colnames(data)[1],"->",colnames(data)[2])
  } else if (F2 > F1) {
    p.value <- 1 - F2/times
    causal <- paste0(colnames(data)[2],"->",colnames(data)[1])
  } else {
    p.value <- NA
    causal <- "not known"
  }
  return(list(f.mean = f.mean, causal.direction = causal, p.value = p.value))
}
