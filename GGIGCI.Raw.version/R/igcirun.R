source("~/Desktop/data_from_HG/IGCI/igci.R")

pair0002 <- read.table("pairs.data/pair0001.txt")
colnames(pair0002) <- c("altitude", "precipitation (yearly value averaged over 1961-1990)")

igci(pair0002$altitude, pair0002$`precipitation (yearly value averaged over 1961-1990)`, refMeasure = "u", estimator = "e")

IGCI.res <- data.frame(pair.name = rep("pair", 108), f1 = rep(0, 108), 
                       f2 = rep(0,108), causal = rep(" ", 108), p = rep(0, 108))
n <- 1
for (i in dir("./pairs.data/", pattern = "[0-9].txt")) {
  message(i)
  pair.name <- strsplit(i, split = ".txt")[[1]][1]
  data.pair <- read.table(paste0("./pairs.data/", i))
  if (dim(data.pair)[2]!=2 || dim(data.pair)[1] < 100) {
    IGCI.res[n,] <- c(pair.name, NA, NA, NA, NA)
  }else {
    #f <- igci(data.pair$V1, data.pair$V2)
    R <- IGCIfor2value(data.pair)
    IGCI.res[n, ] <- c(pair.name, R$f.mean[1], R$f.mean[2], R$causal.direction, R$p.value)
  }
  n <- n+1
  message(paste0(i, " job done!\n"))
}
message(n)

n <- 1
igci.res <- data.frame(pair.name = rep("pair", 108), fue = rep(0, 108), fui = rep(0, 108), fge=rep(0, 108), fgi=rep(0,108), causal = rep("->", 108))
for(i in dir("./pairs.data/", pattern = "[0-9].txt")) {
  message(i)
  pair.name <- strsplit(i, split = ".txt")[[1]][1]
  data.pair <- read.table(paste0("./pairs.data/", i))
  if (dim(data.pair)[2] != 2 || dim(data.pair)[1] < 20) {
    igci.res[n, ] <- c(pair.name, NA, NA)
  } else {
    R <- igci(data.pair[,1],data.pair[,2], refMeasure = "u", estimator = "e")
    igci.res[n,] <- c(pair.name,
                      igci(data.pair$V1, data.pair$V2, refMeasure = "u", estimator = "e"), 
                      igci(data.pair$V1, data.pair$V2, refMeasure = "u", estimator = "i"),
                      igci(data.pair$V1, data.pair$V2, refMeasure = "g", estimator = "e"),
                      igci(data.pair$V1, data.pair$V2, refMeasure = "g", estimator = "i"), "->")
    if (R > 0) {
      igci.res[n, ]$causal <- "<-"
    }
  }
  n <- n+1
  message(paste0(i, " job Done!"))
}

igci.res$causal.r <- pairmeta$causal
write.table(igci.res, file = "igciRresult.txt", quote = F, col.names = T, row.names = F, sep = "\t")

require(devtools)
#install_github("hadley/ggplot2")
#install_github("sachsmc/plotROC")
library(pROC)
data.roc <- na.omit(igci.res)
#data.roc$f <- as.numeric(data.roc$f)
roc.objue <- plot.roc(data.roc$causal.r, as.numeric(data.roc$fue), print.thres = T)
roc.objui <- lines.roc(data.roc$causal.r, as.numeric(data.roc$fui), col = "red")
roc.objge <- lines.roc(data.roc$causal.r, as.numeric(data.roc$fge), col = "green")
roc.objgi <- lines.roc(data.roc$causal.r, as.numeric(data.roc$fgi), col = "blue")
pairmeta <- read.table("pairs.data/pairmeta.txt")
pairmeta$causal <- rep(-1, 108)
pairmeta[which(pairmeta$V2+pairmeta$V3 > pairmeta$V4+pairmeta$V5), ]$causal <- 1
igci.res$causal <- rep("-", 108)
#dim(igci.res)
igci.res$real.f <- rep(0, 108)
igci.res[igci.res$causal == "-",]$real.f <- -1
igci.res[igci.res$causal == "+",]$real.f <- 1
igci.res$f <- as.numeric(igci.res$f)
dim(igci.res51[which((igci.res51$f)*(igci.res51$real.f)>0), ])
igci.res51 <- igci.res[1:51, ]

require(ggplot2)
data.pair <- read.table("./pairs.data/pair0003.txt")
ggplot(data.pair, aes(V1, V2)) +
  geom_point()

