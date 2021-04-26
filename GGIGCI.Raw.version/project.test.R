#setwd("~/Desktop/data_from_HG/my_own_package/GGIGCI/")

source("R/from.assay.to.ped.R")
source("R/read.assay.file.R")
source("R/igci.R")
source("R/read.ped.to.snpMatrix.R")
source("R/impute.res.R")
source("R/GBIGM.without.phen.R")
source("R/GBIGM.R")
source("R/gene.gene.interaction.R")
source("R/gene.gene.causal.R")
require(openssl)

#setwd("../")
data.total.1536 <- read.assay.file(dir = "data_first_20190116/", pattern = "AZ")
selfhead(data.total.1536)
snp.info <- read.xls("~/Desktop/data_from_HG/第一批练习.xls", sheet = 3)
sample.info <- read.xls("~/Desktop/data_from_HG/第一批练习.xls", sheet = 2)

snp.info <- snp.info[which(snp.info$Status=="OK"), ]
sample.info <- sample.info[which(sample.info$Sex!=""), ]
for(i in 1:dim(snp.info)[1]) {
  x <- snp.info$ABI.ID[i]
  snp.info$ABI.ID[i] <- gsub("_+", "_", x)
  #print(x)
}

data.ped <- from.assay.to.ped(data.total.1536, snp.info)
sample.info <- data.frame(Fam.id = sample.info$Unique.Identifier, Sample.Name = sample.info$Unique.Identifier, Paternal.ID = rep(0, dim(sample.info)[1]),
                          Maternal.ID = rep(0, dim(sample.info)[1]), Sex = sample.info$Sex, Phen = sample.info$Classify.)
sample.info[which(sample.info$Sex=="M"), ]$Sex <- 1
sample.info[which(sample.info$Sex=="F"), ]$Sex <- 2
sample.info[which(sample.info$Phen=="NZ"), ]$Phen <- 1
sample.info[which(sample.info$Phen=="SZ"), ]$Phen <- 2

data.ped <- from.assay.to.ped(data.total.1536, snp.info = snp.info, sample.info = sample.info,
                              is.sample.info = T)
write.table(data.ped, file = "data/data.ped", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(data.ped, file = "data/data.with.header.ped", quote = F, sep = "\t", col.names = T, row.names = F)
#res <- read.ped.to.snpMatrix(ped = "data/data.ped", info = snp2gene,
#                             snp.list = snp.list, ped.header = F, ped.sep = "\t")
snp2gene <- read.table("data/data.txt", header = T, sep = "\t")
#The cols of info should be Chr, Gene, Snp and Position.
snp2gene <- snp2gene[,c(3,1,2,4)]
res <- read.ped.to.snpMatrix(ped = "data/data.with.header.ped", ped.header = T, info = snp2gene)

temp.im <- impute.res(res, removed = "sample")
#removed 45 sample in the raw data

#gene2:COMT
#gene1:DLG4

temp.snp.list1 <- temp.im$gene.info[which(temp.im$gene.info$Genenames==gene1), ]$SNPnames
temp.snp.list2 <- temp.im$gene.info[which(temp.im$gene.info$Genenames==gene2), ]$SNPnames
G1 <- temp.im$SnpMatrix[, colnames(temp.im$SnpMatrix)%in%temp.snp.list1]
G2 <- temp.im$SnpMatrix[, colnames(temp.im$SnpMatrix)%in%temp.snp.list2]
phen <- temp.im$sample.info$affected
#GBIGM.without.phen(G1, G2, n.times = 100)
GBIGM(phen, G1, G2, n.times = 100)
temp.ggiud <- gene.gene.interaction(temp.im, n.times = 1000, selected.gene.list = c("DRD1", "DRD2", "DRD3"))
temp.ggc <- gene.gene.causal(temp.ggiud, n.times = 1000, sample.size = 700, fex = F)

table(im1$gene.info$Genenames)
require(igraph)
network.ud <- graph_from_data_frame(temp.ggiud$GGI.ud[which(temp.ggiud$GGI.ud$P.val<0.05), ],
                                    directed = F, vertices = temp.ggiud$GENE.table)
plot(network.ud)

network.d <- graph_from_data_frame(temp.ggc$gg.causal, directed = T, vertices = temp.ggc$gene.table)

plot(network.d)
############test the function of fex #########################
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
    frange <- c(pos.v[1], pos.v[length(pos)])
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

snp.gene1 <- ggiud$gene.info[which(ggiud$gene.info$Genenames==gene1),]
snp.gene2 <- ggiud$gene.info[which(ggiud$gene.info$Genenames==gene2),]
snpm.gene1 <- as(ggiud$SnpMatrix[, which(colnames(ggiud$SnpMatrix)%in%snp.gene1$SNPnames)], "numeric")
snpm.gene2 <- as(ggiud$SnpMatrix[, which(colnames(ggiud$SnpMatrix)%in%snp.gene2$SNPnames)], "numeric")
snpm.pos1 <- snp.gene1$Position
names(snpm.pos1) <- snp.gene1$SNPnames
snpm.pos2 <- snp.gene2$Position
names(snpm.pos2) <- snp.gene2$SNPnames
pos <- snpm.pos1
idx <- order(pos)
pos <- pos[idx]
pos <- (pos - pos[1])/(pos[length(pos)] - pos[1])

fex.gene1 <- fex(snpm.gene1, snp.gene1)
fex.gene2 <- fex(snpm.gene2, snp.gene2)

for(i in 1:3) {
  print(igci(x = fex.gene1[, i], y = fex.gene2[, 1]))
}

for(i in seq_len(dim(selected.gene.pairs)[1])){
  gene1 <- selected.gene.pairs[i,1]
  gene2 <- selected.gene.pairs[i, 2]
  snp.gene1 <- ggiud$gene.info[which(ggiud$gene.info[,2]==gene1),]
  snp.gene2 <- ggiud$gene.info[which(ggiud$gene.info[,2]==gene2),]
  snpm.gene1 <- as(ggiud$SnpMatrix[, colnames(ggiud$SnpMatrix)%in%snp.gene1[,3]], "numeric")
  snpm.gene2 <- as(ggiud$SnpMatrix[, colnames(ggiud$SnpMatrix)%in%snp.gene2[,3]], "numeric")

  #print(dim(fex(snpm.gene1, snp.gene1)))
  #print(dim(fex(snpm.gene2, snp.gene2)))
  fex.gene1 <- fex(snpm.gene1, snp.gene1)
  fex.gene2 <- fex(snpm.gene2, snp.gene2)

  F.matrix <- diag(0, nrow = dim(fex.gene1)[2], ncol = dim(fex.gene2)[2])
  rownames(F.matrix) <- colnames(fex.gene1)
  colnames(F.matrix) <- colnames(fex.gene2)

  for(n in seq_len(ncol(fex.gene1))){
    for(m in seq_len(ncol(fex.gene2))) {
      F.matrix[n, m] <- igci(fex.gene1[, n], fex.gene2[, m])
    }
  }
  print(F.matrix)
}
#########################################
source("./R/gene.gene.causal.R")
gene.gene.causal(ggiud)
#ggc
#ggiud.all
network.ud.all <- graph_from_data_frame(ggiud.all$GGI.ud[ggiud.all$GGI.ud$P.val<=0.1, ], directed = F, vertices = ggiud.all$GENE.table)
network.ud.from.string <- graph_from_data_frame(smaller.network.distance, directed = F, vertices = ggiud.all$GENE.table)
smaller.network.distance$myPvalue <- 0
for(i in seq_len(dim(smaller.network.distance)[1])) {
  #print(my.p[smaller.network.distance$V.from[i], smaller.network.distance$V.to[i]])
  smaller.network.distance$myPvalue[i] <- my.p[smaller.network.distance$V.from[i], smaller.network.distance$V.to[i]]
}
smaller.network.distance$interaction <- 0
smaller.network.distance$interaction <- ifelse(smaller.network.distance$myPvalue <= 0.1, 1, 0)
require(pROC)
roc.my.p.value <- roc(smaller.network.distance$interaction, smaller.network.distance$distance)
plot.roc(roc.my.p.value, print.thres = T)
roc.my.p.value$auc
gbigm.p <- read.table(file = "~/Desktop/GBIGM.p.txt", header = T, sep = "\t")
smaller.network.distance$gbigm.p <- 0
for(i in seq_len(dim(smaller.network.distance)[1])){
  if (!is.null(gbigm.p[smaller.network.distance$V.from[i], smaller.network.distance$V.to[i]])){
    smaller.network.distance$gbigm.p[i] <- gbigm.p[smaller.network.distance$V.from[i], smaller.network.distance$V.to[i]]
  }
  else {
    smaller.network.distance$gbigm.p[i] <- NA
  }
}
smaller.network.distance$interaction.gbigm <- ifelse(is.na(smaller.network.distance$gbigm.p),NA,
                                                     ifelse(smaller.network.distance$gbigm.p < 0.05, 1, 0))
roc.gbigm.p <- roc(smaller.network.distance$interaction.gbigm, 
                   smaller.network.distance$distance)
roc.my.p.value.naomit <- roc(smaller.network.distance.nona$interaction, 
                             smaller.network.distance.nona$distance)

###########################test GBIGM.gene.level.R##############################
#load the data set temp.im include all data
load("C:/Users/liuliangjie19/Desktop/data_from_HG/my_own_package/home.ggigci/ggigci/data/res.after.impute.rda")
#load the gene gene interaction from the GBIGM method
load("C:/Users/liuliangjie19/Desktop/ggiud.data17.rda")
gene.info <- temp.im$gene.info
snp.info <- temp.im$SnpMatrix
#gene1 <- unique(temp.im$gene.info$Genenames)[1]
#gene2 <- unique(temp.im$gene.info$Genenames)[2]
gene1 <- "COMT"
gene2 <- "DLG4"
phen <- temp.im$sample.info$affected
snp.gene1 <- gene.info[which(gene.info[, 2]==gene1), 3]
snp.gene2 <- gene.info[which(gene.info[, 2]==gene2), 3]
snp.info.gene1 <- snp.info[, colnames(snp.info)%in%snp.gene1]
snp.info.gene2 <- snp.info[, colnames(snp.info)%in%snp.gene2]

gguid.gene.level <- gene.gene.interaction(temp.im, gene.level = T)
test.result <- GBIGM.gene.level(phen = phen, gene1 = snp.info.gene1, gene2 = snp.info.gene2, n.times = 1000)
#x.total <- as(snp.info, "numeric")
#apply(x.total[,1:10], 2, table)
x1 <- as(snp.info.gene1, "numeric")
x2 <- as(snp.info.gene2, "numeric")
###the convert method in gbigm.gene.level##
xt1 <- apply(x1, 2, correct.genotype)
xt2 <- apply(x2, 2, correct.genotype)
####
####the convert method in gbigm.normal###
xt1 <- apply(x1, 1, paste, collapse = " ")
xt2 <- apply(x2, 1, paste, collapse = " ")
xt1.2 <- apply(cbind(x1, x2), 1, paste, collapse = " ")
#table(apply(xt1, 1, is.mutation))
#table(apply(xt2, 1, is.mutation))
num.sample <- length(phen)
g1 <- apply(xt1, 1, is.mutation)
g2 <- apply(xt2, 1, is.mutation)
#g1.2 <- apply(cbind(xt1, xt2), 1, is.mutation)
g1.2 <- apply(cbind(g1, g2), 1, paste, collapse = " ")

gx1 <- apply(cbind(phen, g1), 1, paste, collapse = " ")
gx2 <- apply(cbind(phen, g2), 1, paste, collapse = " ")
gx1.2 <- apply(cbind(phen, g1.2), 1, paste, collapse = " ")

for(gene in unique(gene.info$Genenames)){
  snp.gene <- gene.info[gene.info$Genenames==gene, 3]
  snp.info.gene <- snp.info[, colnames(snp.info)%in%snp.gene]
  x <- as(snp.info.gene, "numeric")
  xt <- apply(x, 2, correct.genotype)
  print(gene)
  print(table(apply(xt, 1, is.mutation)))
}
#length(apply(xt1, 1, is.mutation))

gguid.normal$S.matrix
gguid.normal$P.matrix
##############
###########the path of R scripts in the windows system########
#source("C:/Users/liuliangjie19/Desktop/data_from_HG/my_own_package/GGIGCIraw/R/GBIGM.gene.level.R")
#source("C:/Users/liuliangjie19/Desktop/data_from_HG/my_own_package/GGIGCIraw/R/gene.gene.interaction.R")
########compare the gguid result with the netwrok from string database and Hi-c data######

#required data###
V.list.first.comb.h
V.list.first.comb.I
p.matrix.gene.level <- gguid.gene.level$P.matrix
colnames(p.matrix.gene.level)[6] <- "DAOA"
rownames(p.matrix.gene.level)[6] <- "DAOA"
s.matrix.gene.level <- gguid.gene.level$S.matrix
colnames(s.matrix.gene.level)[6] <- "DAOA"
rownames(s.matrix.gene.level)[6] <- "DAOA"
p.matrix.normal <- gguid.normal$P.matrix
colnames(p.matrix.normal)[6] <- "DAOA"
rownames(p.matrix.normal)[6] <- "DAOA"
s.matrix.normal <- gguid.normal$S.matrix
colnames(s.matrix.normal)[6] <- "DAOA"
rownames(s.matrix.normal)[6] <- "DAOA"
#plot the 17 gene in 3D####
library(scatterplot3d)
library(rgl)
coord.h.17 <- coord.h[coord.h$gene.symbol%in%colnames(p.matrix.normal), ]
coord.h.17$color <- ifelse(coord.h.17$Chr=="X", 24, coord.h.17$Chr)
coord.I.17 <- coord.I[coord.I$gene.symbol%in%colnames(p.matrix.normal), ]
coord.I.17$color <- ifelse(coord.I.17$Chr=="X", 24, coord.I.17$Chr)
#palette.colors(n = 25)
scatterplot3d(coord.h.17[, 6:8], color = coord.h.17$color, type = 'p',highlight.3d=F,
              angle=60,grid=T,box=T,scale.y=1, cex.symbols=1.2, pch=16)
scatterplot3d(coord.I.17[, 6:8], color = coord.I.17$color, type = 'p',highlight.3d=F,
              angle=60,grid=T,box=T,scale.y=1, cex.symbols=1.2, pch=16)
#palette("default")
d.list <- NULL
for(i in 1:17){
  d <- sqrt(sum((coord.h.17[i, 6:8]-coord.I.17[i, 6:8])^2))
  d.list <- c(d.list, d)
}
boxplot(d.list)
#################
compare.gguid.network <- function(v.list, p.matrix, s.matrix){
  compare.result <- v.list
  compare.result$p.value <- NA
  compare.result$s.value <- NA
  if (!all(dim(p.matrix)==dim(s.matrix))){
    stop("p.matrix and s.matrix must be from the same GGUID result.")
  }
  if (!all(unique(compare.result$gene1)%in%colnames(p.matrix))){
    stop("the V.list and the matrixes must be from the same gene list.")
  }
  n <- dim(v.list)[1]
  bar <- txtProgressBar(0, n, char = ">", style = 3)
  for(i in seq_len(n)){
    from <- compare.result[i, 1]
    to <- compare.result[i, 2]
    compare.result$p.value[i] <- p.matrix[from, to]
    compare.result$s.value[i] <- s.matrix[from, to]
    
    setTxtProgressBar(bar, i)
  }
  
  return(compare.result)
}


compare.gene.level.to.h <- compare.gguid.network(V.list.first.comb.h, 
                                                 p.matrix.gene.level, s.matrix.gene.level)
compare.gene.level.to.h[,4] <- as(compare.gene.level.to.h[, 4], "numeric")
summary(compare.gene.level.to.h)
compare.gene.level.to.I <- compare.gguid.network(V.list.first.comb.I, 
                                                 p.matrix.gene.level, s.matrix.gene.level)
compare.gene.level.to.I[, 4] <- as(compare.gene.level.to.I[, 4], "numeric")
summary(compare.gene.level.to.I)


compare.normal.to.h <- compare.gguid.network(V.list.first.comb.h, 
                                             p.matrix.normal, s.matrix.normal)
compare.normal.to.h[, 4] <- as(compare.normal.to.h[, 4], "numeric")
summary(compare.normal.to.h)
compare.normal.to.I <- compare.gguid.network(V.list.first.comb.I, 
                                             p.matrix.normal, s.matrix.normal)
compare.normal.to.I[, 4] <- as(compare.normal.to.I[, 4], "numeric")
summary(compare.normal.to.I)
compare.normal.to.h$uid <- ifelse(compare.normal.to.h$p.value<0.05, 1, 0)
compare.normal.to.I$uid <- ifelse(compare.normal.to.I$p.value<0.05, 1, 0)
compare.gene.level.to.h$uid <- ifelse(compare.gene.level.to.h$p.value<0.05, 1, 0)
compare.gene.level.to.I$uid <- ifelse(compare.gene.level.to.I$p.value<0.05, 1, 0)
#the roc result of the total 17 genes ####
roc.normal.to.h <- roc(compare.normal.to.h$uid, compare.normal.to.h$distance.weight)
#plot the roc plot between GBIGM result and distance#
get.auc.table <- function(gene.table, compare.h, compare.I, k = 1){
  auc.list.weight <- NULL
  auc.list.coordinate.h <- NULL
  auc.list.coordinate.I <- NULL
  for (i in k:8){
    gene.most.loci <- gene.table[gene.table$Snps.Freq>i, 1]
    #correct the gene name
    gene.most.loci <- ifelse(gene.most.loci=="G72/G30", "DAOA", gene.most.loci)
    temp.compare <- compare.h[compare.h[,1]%in%gene.most.loci & compare.h[,2]%in%gene.most.loci, ]
    roc.temp.weight <- roc(temp.compare$uid, temp.compare$distance.weight)
    #roc.temp.weight <- roc(temp.compare$uid, temp.compare$distance.weight, levels = c("1", "0"), direction = "<")
    roc.temp.coordinate <- roc(temp.compare$uid, temp.compare$distance.coordinate)
    #roc.temp.coordinate <- roc(temp.compare$uid, temp.compare$distance.coordinate, levels = c("1", "0"), direction = "<")
    auc.list.weight <- c(auc.list.weight, auc(roc.temp.weight))
    auc.list.coordinate.h <- c(auc.list.coordinate.h, auc(roc.temp.coordinate))
    temp.compare <- compare.I[compare.I[,1]%in%gene.most.loci & compare.I[,2]%in%gene.most.loci, ]
    roc.temp.coordinate <- roc(temp.compare$uid, temp.compare$distance.coordinate)
    #roc.temp.coordinate <- roc(temp.compare$uid, temp.compare$distance.coordinate, levels = c("1", "0"), direction = "<")
    auc.list.coordinate.I <- c(auc.list.coordinate.I, auc(roc.temp.coordinate))
  }
  auc.temp <- data.frame(n = c(1:8), 
                         coordinate.h = auc.list.coordinate.h, 
                         coordinate.I = auc.list.coordinate.I, 
                         weight = auc.list.weight)
  return(auc.temp)
}
#when i==6 the auc > 0.9
auc.temp.normal <- get.auc.table(gene.table = gguid.normal$GENE.table, 
                                 compare.h = compare.normal.to.h,
                                 compare.I = compare.normal.to.I)
auc.temp.gene.level <- get.auc.table(gene.table = gguid.gene.level$GENE.table, 
                                     compare.h = compare.gene.level.to.h,
                                     compare.I = compare.gene.level.to.I)
plot(auc.temp$n, auc.temp$weight, type = "l", ylim = c(0.4, 1), xlab = "number of snp", ylab = "AUC")
lines(auc.temp$n, auc.temp$coordinate.I, col = "red")
lines(auc.temp$n, auc.temp$coordinate.h, col = "blue")
ggplot(data = auc.temp) + geom_line(aes(x = n, y = weight, color = "weight")) + 
  geom_line(aes(x = n, y = coordinate.I, color = "coordinate.IMR")) +
  geom_line(aes(x = n, y = coordinate.h, color = "coordinate.hESC")) + 
  theme_classic() + ylab("AUC") + ylim(0.4, 1)


#plot(roc.temp.weight, print.thres = T)
#the temp.compare is under the i==6###
p.heatmap.temp.compare <- ggplot(temp.compare,  aes(gene1, gene2)) +
  geom_tile(aes(fill = p.value), colour = "white")
p.heatmap.temp.compare + scale_fill_gradient(name = "p.value", low = "red", high = "white") +
  geom_text(aes(gene1, gene2, label = p.value))
roc.i.6 <- roc(temp.compare$uid, temp.compare$distance.weight)
plot(roc.i.6, print.thres = T)
############use ggplot2 to plot the heatmap for the gguid####################
#######the data required######
gguid.gene.level
gguid.normal
#########the example code for plotting the heatmap#####################
require(reshape2)
p.heatmap.gene.level <- ggplot(gguid.gene.level$GGI.ud, aes(from, to)) +
  geom_tile(aes(fill = P.val), colour = "white")
p.heatmap.gene.level + scale_fill_gradient(name="P.value", low = "red",high = "white") +
  geom_text(aes(from, to, label = P.val), color = "black", size = 3)

p.heatmap.normal <- ggplot(gguid.normal$GGI.ud, aes(from, to)) +
  geom_tile(aes(fill = P.val), colour = "white")
p.heatmap.normal + scale_fill_gradient(name = "P.value", low = "red", high = "white") +
  geom_text(aes(from, to, label = P.val), color = "black", size = 3)

p.heatmap.normal.multi.loci <- ggplot(gguid.normal.multi.loci$GGI.ud, aes(from, to)) +
  geom_tile(aes(fill = P.val), colour = "white")
p.heatmap.normal.multi.loci + scale_fill_gradient(name = "P.value", low = "red", high = "white") +
  geom_text(aes(from, to, label = P.val), color = "black", size = 3)
##############################################################################

#########filter the gene with only one or two snp loci#########################
gene.multi.loci <- gguid.normal$GENE.table[gguid.normal$GENE.table$Snps.Freq>6, 1]
#run as the local job ###
gguid.gene.level.multi.loci <- gene.gene.interaction(temp.im, gene.level = T, selected.gene.list = gene.multi.loci)
gguid.normal.multi.loci <- gene.gene.interaction(temp.im, gene.level = F, selected.gene.list = gene.multi.loci)
#########################

######the gguid result by using the other methods################################
source("~/Desktop/GGIGCI.Raw.version/GGI.other.methods/selectSnps.R")
source("~/Desktop/GGIGCI.Raw.version/GGI.other.methods/GGI.R")
source("~/Desktop/GGIGCI.Raw.version/GGI.other.methods/minP.R")
source("~/Desktop/GGIGCI.Raw.version/GGI.other.methods/CCA.R")
source("~/Desktop/GGIGCI.Raw.version/GGI.other.methods/KCCA.R")
source("~/Desktop/GGIGCI.Raw.version/GGI.other.methods/CLD.R")
source("~/Desktop/GGIGCI.Raw.version/GGI.other.methods/PLSPM.R")
source("~/Desktop/GGIGCI.Raw.version/GGI.other.methods/GATES.R")
source("~/Desktop/GGIGCI.Raw.version/GGI.other.methods/tTS.R")
source("~/Desktop/GGIGCI.Raw.version/GGI.other.methods/tProd.R")
source("~/Desktop/GGIGCI.Raw.version/GGI.other.methods/PCA.R")
source("~/Desktop/GGIGCI.Raw.version/GGI.other.methods/SSIBase.R")
#require data temp.im
gene.most.loci <- gguid.normal$GENE.table[gguid.normal$GENE.table$Snps.Freq>2, 1]
temp.im.other.methods <- selectSnps(temp.im$SnpMatrix, temp.im$gene.info, select = gene.most.loci)
gguid.normal.minp <- GGI(phen = phen, snpmatrix = temp.im.other.methods$snpX, 
                        genes.info = temp.im.other.methods$genes.info, 
                        method = "minP")
gguid.normal.CCA <- GGI(phen = phen, snpmatrix = temp.im.other.methods$snpX, 
                        genes.info = temp.im.other.methods$genes.info,
                        method = "CCA")
gguid.normal.KCCA <- GGI(phen = phen, snpmatrix = temp.im.other.methods$snpX,
                         genes.info = temp.im.other.methods$genes.info, 
                         method = "KCCA")
gguid.normal.CLD <- GGI(phen = phen, snpmatrix = temp.im.other.methods$snpX,
                        genes.info = temp.im.other.methods$genes.info,
                        method = "CLD")
gguid.normal.PLSPM <- GGI(phen = phen, snpmatrix = temp.im.other.methods$snpX,
                          genes.info = temp.im.other.methods$genes.info,
                          method = "PLSPM")
gguid.normal.GATES <- GGI(phen = phen, snpmatrix = temp.im.other.methods$snpX,
                          genes.info = temp.im.other.methods$genes.info,
                          method = "GATES")
gguid.normal.tTS <- GGI(phen = phen, snpmatrix = temp.im.other.methods$snpX,
                        genes.info = temp.im.other.methods$genes.info,
                        method = "tTS")
gguid.normal.tProd <- GGI(phen = phen, snpmatrix = temp.im.other.methods$snpX,
                          genes.info = temp.im.other.methods$genes.info,
                          method = "tProd")
####################################################################################
###CCA minp plspm tprod tts###
#cca
cca.result <- melt(gguid.normal.CCA$p.value)
colnames(cca.result) <- c("gene1", "gene2", "cca.p")
minp.result <- melt(gguid.normal.minp$p.value)
colnames(minp.result) <- c("gene1", "gene2", "minp.p")
#plspm.result <- melt(gguid.normal.PLSPM$p.value)
#colnames(plspm.result) <- c("gene1", "gene2", "plspm.p")
tprod.result <- melt(gguid.normal.tProd$p.value)
colnames(tprod.result) <- c("gene1", "gene2", "tprod.p")
tts.result <- melt(gguid.normal.tTS$p.value)
colnames(tts.result) <- c("gene1", "gene2", "tts.p")
compare.other.method.to.h <- na.omit(merge(compare.normal.to.h, cca.result, by = c("gene1", "gene2"), all.x = T))
compare.other.method.to.h$uid.cca <- ifelse(compare.other.method.to.h$cca.p<0.05, 1, 0)
compare.other.method.to.h <- na.omit(merge(compare.other.method.to.h, minp.result, by = c("gene1", "gene2"), all.x = T))
compare.other.method.to.h$uid.minp <- ifelse(compare.other.method.to.h$minp.p<0.05, 1, 0)
#compare.other.method.to.h <- na.omit(merge(compare.other.method.to.h, plspm.result, by = c("gene1", "gene2"), all.x = T))
#compare.other.method.to.h$uid.plspm <- ifelse(compare.other.method.to.h$plspm.p<0.05, 1, 0)
compare.other.method.to.h <- na.omit(merge(compare.other.method.to.h, tprod.result, by = c("gene1", "gene2"), all.x = T))
compare.other.method.to.h$uid.tprod <- ifelse(compare.other.method.to.h$tprod.p<0.05, 1, 0)
compare.other.method.to.h <- na.omit(merge(compare.other.method.to.h, tts.result, by = c("gene1", "gene2"), all.x = T))
compare.other.method.to.h$uid.tts <- ifelse(compare.other.method.to.h$tts.p<0.05, 1, 0)
get.auc.table.other.method <- function(gene.table, compare.other.method, weight = 3, k.min = 2, k.max = 8){
  auc.list <- NULL
  compare.other.method <- compare.other.method[, c(1,2,weight, grep("uid", colnames(compare.other.method)))]
  for (i in k.min:k.max) {
    gene.most.loci <- gene.table[gene.table[,2]>i, 1]
    temp.compare <- compare.other.method[compare.other.method[,1]%in%gene.most.loci & compare.other.method[,2]%in%gene.most.loci, ]
    for (j in 4:dim(temp.compare)[2]){
      #roc.temp.weight <- roc(as.factor(temp.compare[,j]), temp.compare[,3], levels = c("1", "0"), direction = "<")
      roc.temp.weight <- roc(temp.compare[,j], temp.compare[, 3])
      auc.list <- c(auc.list, auc(roc.temp.weight))
    }
  }
  l <- k.max - k.min + 1
  return(matrix(auc.list, nrow = l, byrow = T))
}
auc.table.compare.with.other.method <- get.auc.table.other.method(gguid.normal$GENE.table, compare.other.method.to.h, k.max = 5)
auc.table.compare.with.other.method <- get.auc.table.other.method(gguid.normal$GENE.table, compare.other.method.to.h, weight = 4, k.max = 5)
colnames(auc.table.compare.with.other.method) <- c("GBIGM", "cca", "minp", "tprod", "tts")
auc.table.compare.with.other.method <- as.data.frame(auc.table.compare.with.other.method)
auc.table.compare.with.other.method$n <- 2:5
ggplot(data = auc.table.compare.with.other.method) + geom_line(aes(x = n, y = GBIGM, color = "GBIGM")) +
  geom_line(aes(x = n, y = cca, color = "cca")) +
  geom_line(aes(x = n, y = minp, color = "minp")) +
  geom_line(aes(x = n, y = tprod, color = "tprod")) + 
  geom_line(aes(x = n, y = tts, color = "tts")) + theme_classic() + ylab("AUC")
#compare.other.method.to.I <- na.omit(merge(compare.normal.to.I, cca.result, by = c("gene1", "gene2"), all.x = T))

###########plot the network graph for the interaction from string database###########
#require data-set: interaction.136, gene.info and all.protein###
summary(interaction.136)
genename.17 <- unique(gene.info$Genenames)
genename.17[6] <- "DAOA"
temp <- all.protein[all.protein$gene.symbol%in%genename.17, 1:2]
genename.17 <- temp$gene.symbol
names(genename.17) <- temp$ensp.id
rm(temp)
gene1 <- genename.17[interaction.136$protein1]
gene2 <- genename.17[interaction.136$protein2]
interaction.136$gene1 <- gene1
interaction.136$gene2 <- gene2
interaction.136 <- interaction.136[,c(17,18, 1:16)]

require(igraph)
require(ggraph)
require(tidygraph)
require(tidyverse)
v.table <- gguid.normal$GENE.table
v.table[6,1] <- "DAOA"
for(i in 5:18) {
  if (all(interaction.136[, i]==0)){
    next
  }
  plot.title <- colnames(interaction.136)[i]
  temp <- interaction.136[, c(1,2,i)]
  colnames(temp)[3] <- "weight"
  #temp.v <- unique(c(temp[,1],temp[,2]))
  temp <- temp[temp[,3]!=0, ]
  temp.graph <- graph_from_data_frame(temp, directed = F, vertices = v.table)
  #plot(temp.graph)
  temp.graph <- as_tbl_graph(temp.graph) %>% mutate(deg = centrality_degree())
  print(plot.title)
  #autograph(temp.graph)
  #pdf(file = paste0("~/Desktop/data_from_HG/string/first_data_17_gene/", plot.title, ".pdf"))
  ggraph(temp.graph, layout = "kk") + 
    geom_edge_fan(aes(edge_width=weight),color="lightblue", show.legend = T, end_cap = circle(1, 'mm')) +
    geom_node_point(aes(size = deg), fill = "green", shape=21) +
    geom_node_text(aes(label = name), size= 2.5, repel = T) +
    scale_edge_width_continuous(range = c(0.2,1)) +
    theme_graph() + theme(text = element_text(family = "mono"))
  ggsave(file = paste0("~/Desktop/data_from_HG/string/first_data_17_gene/", plot.title, ".pdf"))
  #print(temp.plot)
  #dev.off()
  #break
}
rm(temp)
##############################################################################################

###########plot the network of my 17 genes network from GBIGM using ggraph########
temp.e <- gguid.normal$GGI.ud
temp.v <- gguid.normal$GENE.table
temp.v[,2] <- as.numeric(temp.v[,2])
temp.e[which(temp.e$from=="G72/G30"), 1] <- "DAOA"
temp.e[which(temp.e$to=="G72/G30"), 2] <- "DAOA"
temp.v[which(temp.v$Genenames=="G72/G30"), 1] <- "DAOA"
#for(i in 3:4) {
  i <- 3
  plot.title <- colnames(temp.e)[i]
  temp <- temp.e[, c(1,2,i)]
  if (i==3){
    temp[, 3] <- ifelse(temp[,3]==0, -log(0.0001), -log(temp[,3]))
    temp <- temp[which(temp[, 3]>=-log(0.05)), ]
  }
  colnames(temp)[3] <- "weight"
  temp.graph <- graph_from_data_frame(temp, directed = F, vertices = temp.v)
  temp.graph <- as_tbl_graph(temp.graph) %>% mutate(deg = centrality_degree())
  print(plot.title)
  temp.layout <- layout_with_fr(temp.graph)
  ggraph(temp.graph, layout = temp.layout) +
    geom_edge_fan(aes(edge_width=weight),color="lightblue", show.legend = T, end_cap = circle(1, 'mm')) +
    geom_node_point(aes(size = Snps.Freq, fill = 11 - deg), shape=21) +
    geom_node_text(aes(label = name), size= 2.5, repel = T) +
    scale_edge_width_continuous(range = c(0.2,1)) +
    theme_graph() + theme(text = element_text(family = "mono"))
  ggsave(file = paste0("~/Desktop/data_from_HG/string/first_data_17_gene/", plot.title, ".pdf"))
#}
#############################################################################################
##calculate the gene gene causal betweeen 17 gene #require data: gguid.normal gguid.normal.multi.loci
ggc.all.17 <- gene.gene.causal(data.ggiud = gguid.normal, fex = F)
ggc.multi.loci <- gene.gene.causal(data.ggiud = gguid.normal.multi.loci, fex = F)
ggc.all.17.fex <- gene.gene.causal(data.ggiud = gguid.normal)
ggc.multi.loci.fex <- gene.gene.causal(data.ggiud = gguid.normal.multi.loci)
##use ggraph to plot ggc network between 17 genes or genes with multi locis.##
#ggc.all
E.temp <- ggc.all.17$gg.causal
E.temp$from[which(E.temp$from=="G72/G30")] <- "DAOA"
V.temp <- ggc.all.17$gene.table
V.temp$Genenames[which(V.temp$Genenames=="G72/G30")] <- "DAOA"

E.temp <- E.temp[which(ggc.all.17$gg.f$P.value<=0.05), ]
#E.temp$F <- ifelse(E.temp$causal.direction=="gene1->gene2", abs(as.numeric(E.temp$f1)), abs(as.numeric(E.temp$f2)))
E.temp$F <- abs(as.numeric(E.temp$F.value))
require(igraph)
require(ggraph)
temp.graph <- graph_from_data_frame(d = E.temp, directed = T, vertices = V.temp)
temp.graph <- as_tbl_graph(temp.graph) %>% mutate(deg.out = centrality_degree(mode = "out"), deg.in = centrality_degree(mode = "in"))
layout.temp <- layout_with_fr(temp.graph)
ggraph(temp.graph, layout = "centrality", cent = deg.out)+
  geom_edge_fan(aes(start_cap = label_rect(from), end_cap = label_rect(to), edge_width = F), 
                color = "lightblue", show.legend = T, arrow = arrow(length=unit(2, 'mm'))) +
  geom_node_text(aes(label = name), size = 6) +
  scale_edge_width_continuous(range = c(0.2, 1)) + 
  theme_graph() + theme(text = element_text(family = "mono"))

########
#ggc.multi.loci
E.temp <- ggc.multi.loci$gg.causal
#E.temp$from[which(E.temp$from=="G72/G30")] <- "DAOA"
V.temp <- ggc.multi.loci$gene.table
#V.temp$Genenames[which(V.temp$Genenames=="G72/G30")] <- "DAOA"

E.temp <- E.temp[which(ggc.multi.loci$gg.f$P.value<=0.05), ]
#E.temp$F <- ifelse(E.temp$causal.direction=="gene1->gene2", abs(as.numeric(E.temp$f1)), abs(as.numeric(E.temp$f2)))
E.temp$F <- abs(as.numeric(E.temp$F.value))
require(igraph)
require(ggraph)
temp.graph <- graph_from_data_frame(d = E.temp, directed = T, vertices = V.temp)
temp.graph <- as_tbl_graph(temp.graph) %>% mutate(deg.out = centrality_degree(mode = "out"), deg.in =centrality_degree(mode = "in"))
layout.temp <- layout_with_fr(temp.graph)
ggraph(temp.graph, layout = layout.temp)+
  geom_edge_fan(aes(start_cap = label_rect(from), end_cap = label_rect(to), edge_width = F), 
                color = "lightblue", show.legend = T, arrow = arrow(length=unit(2, 'mm'))) +
  geom_node_text(aes(label = name), size = 6) +
  scale_edge_width_continuous(range = c(0.2, 1)) + 
  theme_graph() + theme(text = element_text(family = "mono"))

###compare the GBIGM result with the data from string and Hi-C and other methods' result###############################
temp <- gguid.normal$GGI.ud
temp.interaction.136 <- interaction.136
temp.interaction.136$p.val <- rep(0, dim(temp.interaction.136)[1])
temp.interaction.136$s.val <- rep(0, dim(temp.interaction.136)[1])

for(i in 1:dim(temp)[1]) {
  gene.pairs <- c(temp[i,1], temp[i, 2])
  p.s <- temp[i, 3:4]
  temp.interaction.136[temp.interaction.136[,1]%in%gene.pairs & temp.interaction.136[, 2]%in%gene.pairs, 19:20] <- p.s
}
interaction.136.pval.sval <- temp.interaction.136
rm(temp.interaction.136)
temp.interaction.136 <- interaction.136.pval.sval[, 5:20]
temp.interaction.136 <- cbind(interaction.136.pval.sval[, 1:4], temp.interaction.136[,colSums(interaction.136.pval.sval[,5:20])!=0])
require(ggpubr)
#ggdotplot(temp.interaction.136, x = "s.val", y = "p.val", size = 1, binwidth = 0.02, add = c("violin", "mean_sd"))
cor(temp.interaction.136[,5:16])
require(corrplot)
corrplot(cor(temp.interaction.136[, 5:16]), order = "AOE", type = "upper", tl.pos = "tp")
corrplot(cor(temp.interaction.136[, 5:16]), add = T, type = "lower", method = "number", order="AOE", col="black",diag=FALSE,tl.pos="n", cl.pos="n")

temp.gene.multi.loci <- temp.interaction.136[temp.interaction.136[, 1]%in%gene.multi.loci & temp.interaction.136[, 2]%in%gene.multi.loci, ]
cor(temp.gene.multi.loci[, 5:16])
corrplot(cor(temp.gene.multi.loci[, 5:16]), order = "AOE", type = "upper", tl.pos = "tp")
corrplot(cor(temp.interaction.136[, 5:16]), add = T, type = "lower", method = "number", order = "AOE", col = "black", diag=F, tl.pos = "n", cl.pos = "n")

gene.distance.136 <- read.table("~/Desktop/gene.distance.tsv")
gene.distance.136 <- gene.distance.136[, 2:5]
colnames(gene.distance.136) <- c("protein1", "protein2", "step", "distance")
pair.136.compare <- gene.distance.136
x <- gguid.normal.default$P.matrix
y <- gguid.normal.default$S.matrix
pair.136.compare$distance.I <- NA
pair.136.compare$distance.H <- NA
pair.136.compare$p.val <- NA
pair.136.compare$s.val <- NA

colnames(x)[6] <- rownames(x)[6] <- "DAOA"
colnames(y)[6] <- rownames(y)[6] <- "DAOA"
V.list.first.comb.h
V.list.first.comb.I
for(i in 1:nrow(gene.distance.136)) {
  gene1 <- gene.distance.136[i,1]
  gene2 <- gene.distance.136[i,2]
  pair.136.compare$p.val[i] <- x[gene1, gene2]
  pair.136.compare$s.val[i] <- y[gene1, gene2]
  pair.136.compare$distance.I[i] <- as.numeric(V.list.first.comb.I$distance.coordinate[which(sapply(1:nrow(V.list.first.comb.I), 
                                                                                                    FUN = function(k){return(all(V.list.first.comb.I[k, 1:2]%in%c(gene1, gene2)))}))])
  pair.136.compare$distance.H[i] <- as.numeric(V.list.first.comb.h$distance.coordinate[which(sapply(1:nrow(V.list.first.comb.h),
                                                                                                    FUN = function(k){return(all(V.list.first.comb.h[k, 1:2]%in%c(gene1, gene2)))}))])
}

get.auc.table.default <- function(gene.table, pair.compare, k = 1){
  auc.list.weight <- NULL
  auc.list.coordinate.h <- NULL
  auc.list.coordinate.I <- NULL
  for (i in k:8){
    gene.most.loci <- gene.table[gene.table$Snps.Freq>=i, 1]
    #correct the gene name
    gene.most.loci <- ifelse(gene.most.loci=="G72/G30", "DAOA", gene.most.loci)
    temp.compare <- pair.compare[pair.compare[,1]%in%gene.most.loci & pair.compare[,2]%in%gene.most.loci, ]
    roc.temp.weight <- roc(temp.compare$uid, temp.compare$distance)
    #roc.temp.weight <- roc(temp.compare$uid, temp.compare$distance.weight, levels = c("1", "0"), direction = "<")
    roc.temp.coordinate <- roc(temp.compare$uid, temp.compare$distance.H)
    #roc.temp.coordinate <- roc(temp.compare$uid, temp.compare$distance.coordinate, levels = c("1", "0"), direction = "<")
    auc.list.weight <- c(auc.list.weight, auc(roc.temp.weight))
    auc.list.coordinate.h <- c(auc.list.coordinate.h, auc(roc.temp.coordinate))
    roc.temp.coordinate <- roc(temp.compare$uid, temp.compare$distance.I)
    #roc.temp.coordinate <- roc(temp.compare$uid, temp.compare$distance.coordinate, levels = c("1", "0"), direction = "<")
    auc.list.coordinate.I <- c(auc.list.coordinate.I, auc(roc.temp.coordinate))
  }
  auc.temp <- data.frame(n = c(1:8), 
                         coordinate.h = auc.list.coordinate.h, 
                         coordinate.I = auc.list.coordinate.I, 
                         weight = auc.list.weight)
  return(auc.temp)
}

get.auc.table.default(gguid.normal.default$GENE.table, pair.136.compare, k= 1)

pair.136.compare$p.adjsut <- p.adjust(pair.136.compare$p.val, method = "BH")
pair.136.compare$uid <- ifelse(pair.136.compare$p.val<=0.05, 1, 0)
p.heatmap.normal.default <- ggplot(gguid.normal.default$GGI.ud, aes(from, to)) +
  geom_tile(aes(fill = P.val), colour = "white")
p.heatmap.normal.default + scale_fill_gradient(name="P.value", low = "red",high = "white") +
  geom_text(aes(from, to, label = P.val), color = "black", size = 3)

pair.136.compare$cca.p.value <- sapply(1:136, FUN = function(i){return(ifelse(all(c(pair.136.compare[i, 1], pair.136.compare[i, 2])%in%colnames(gguid.normal.CCA$p.value)), 
                                              gguid.normal.CCA$p.value[pair.136.compare[i,1], pair.136.compare[i,2]], NA))})
pair.136.compare$minp.p.value <- sapply(1:136, FUN=function(i){return(ifelse(all(c(pair.136.compare[i, 1], pair.136.compare[i,2])%in%colnames(gguid.normal.minp$p.value)), 
                                                                             gguid.normal.minp$p.value[pair.136.compare[i,1], pair.136.compare[i,2]], NA))})
#pair.136.compare$plspm.p.value <- sapply(1:136, FUN=function(i){return(ifelse(all(c(pair.136.compare[i, 1], pair.136.compare[i,2])%in%colnames(gguid.normal.PLSPM$p.value)),
#                                                                              gguid.normal.PLSPM$p.value[pair.136.compare[i,1], pair.136.compare[i,2]], NA))})
pair.136.compare$tpord <- sapply(1:136, FUN = function(i){return(ifelse(all(c(pair.136.compare[i,1],pair.136.compare[i,2])%in%colnames(gguid.normal.tProd$p.value)), 
                                                                        gguid.normal.tProd$p.value[pair.136.compare[i,1], pair.136.compare[i,2]],NA))})
pair.136.compare$tts <- sapply(1:136, FUN = function(i){return(ifelse(all(c(pair.136.compare[i,1],pair.136.compare[i,2])%in%colnames(gguid.normal.tTS$p.value)), 
                                                                      gguid.normal.tTS$p.value[pair.136.compare[i,1], pair.136.compare[i,2]],NA))})
pair.136.compare$distance.shortest <- NA
for(i in 1:136){
  gene1 <- pair.136.compare[i, 1]
  gene2 <- pair.136.compare[i, 2]
  pair.136.compare$distance.shortest[i] <- V.list.first.comb.h$distance.weight[which(sapply(1:136, FUN = function(k){return(all(V.list.first.comb.h[k,1:2]%in%c(gene1, gene2)))}))]
}

pair.136.compare$combine.score <- NA
for(i in 1:136){
  gene1 <- pair.136.compare[i, 1]
  gene2 <- pair.136.compare[i, 2]
  pair.136.compare$combine.score[i] <- interaction.136$combined_score[which(sapply(1:136, FUN=function(k){return(all(interaction.136[k, 1:2]%in%c(gene1, gene2)))}))]
}

temp.auc.array <- 1:225
dim(temp.auc.array) <- c(5, 5, 9)

for(i in 1:9){
  gene.list <- gguid.normal.default$GENE.table$Genenames[which(gguid.normal.default$GENE.table$Snps.Freq>=i)]
  temp.pair <- pair.136.compare[pair.136.compare$protein1%in%gene.list & pair.136.compare$protein2%in%gene.list, ]
  temp.auc <- NULL
  for (d in c(4:6,14,15)){
    temp.distance <- temp.pair[, d]
    for(p in c(7,10:13)){
      temp.p.val <- temp.pair[, p]
      if (any(is.na(temp.p.val))){
        temp.auc <- c(temp.auc, NA)
        next
      }
      temp.uid <- as.factor(ifelse(temp.p.val<=0.05, 1, 0))
      if (nlevels(temp.uid)==1){
        temp.auc <- c(temp.auc, NA)
        next
      }
      temp.auc <- c(temp.auc, auc(roc(temp.uid, temp.distance)))
    }
  }
  print(matrix(temp.auc, nrow = 5, dimnames = list(c(colnames(temp.pair)[c(7,10:13)]), c(colnames(temp.pair)[c(4:6,14,15)]))))
  #table.name <- paste0("auc.table", i)
  temp.auc.array[,,i] <- matrix(temp.auc, nrow = 5, dimnames = list(c(colnames(temp.pair)[c(7,10:13)]), c(colnames(temp.pair)[c(4:6,14,15)])))
}

###calculate the ggc realtionship basing on gguid.normal.default###
ggc.all.17.default.ge <- gene.gene.causal(gguid.normal.default, fex = F, ref = "Gaussian",est = "entropy")
ggc.all.17.default.gi <- gene.gene.causal(gguid.normal.default, fex = F, ref = "Gaussian", est = "intergral")
ggc.all.17.default.ue <- gene.gene.causal(gguid.normal.default, fex = F, ref = "uniform", est = "entropy")
ggc.all.17.default.ui <- gene.gene.causal(gguid.normal.default, fex = F, ref = "uniform", est = "intergral")
ggc.all.17.default.fex.ge <- gene.gene.causal(gguid.normal.default, fex = T, ref = "Gaussian", est = "entropy")
ggc.all.17.default.fex.gi <- gene.gene.causal(gguid.normal.default, fex = T, ref = "Gaussian", est = "intergral")
ggc.all.17.default.fex.ue <- gene.gene.causal(gguid.normal.default, fex = T, ref = "uniform", est = "entropy")
ggc.all.17.default.fex.ui <- gene.gene.causal(gguid.normal.default, fex = T, ref = "uniform", est = "intergral")

compare.between.ggc.evec$snp.freq <- gguid.normal.default$GENE.table$Snps.Freq
gene.list <- compare.between.ggc.evec$gene.symbol
gene.list[6] <- "G72/G30"

get.degree.in.gcc <- function(gene.list, gg.causal, mode = "from"){
  if (mode=="from"){
    return(sapply(gene.list, FUN = function(gene){return(sum(gg.causal$from==gene))}))
  }
  else if (mode=="to") {
    return(sapply(gene.list, FUN = function(gene){return(sum(gg.causal$to==gene))}))
  }
}
get.degree.in.gcc(gene.list, ggc.all.17.default.ge$gg.causal, mode = "from")
temp.degree.ggc <- cbind(get.degree.in.gcc(gene.list, ggc.all.17.default.ge$gg.causal, mode = "from"), 
      get.degree.in.gcc(gene.list, ggc.all.17.default.ge$gg.causal, mode="to"),
      get.degree.in.gcc(gene.list, ggc.all.17.default.gi$gg.causal, mode="from"),
      get.degree.in.gcc(gene.list, ggc.all.17.default.gi$gg.causal, mode="to"),
      get.degree.in.gcc(gene.list, ggc.all.17.default.ue$gg.causal, mode="from"),
      get.degree.in.gcc(gene.list, ggc.all.17.default.ue$gg.causal, mode="to"),
      get.degree.in.gcc(gene.list, ggc.all.17.default.ui$gg.causal, mode="from"),
      get.degree.in.gcc(gene.list, ggc.all.17.default.ui$gg.causal, mode="to"),
      get.degree.in.gcc(gene.list, ggc.all.17.default.fex.ge$gg.causal, mode="from"),
      get.degree.in.gcc(gene.list, ggc.all.17.default.fex.ge$gg.causal, mode="to"),
      get.degree.in.gcc(gene.list, ggc.all.17.default.fex.gi$gg.causal, mode="from"),
      get.degree.in.gcc(gene.list, ggc.all.17.default.fex.gi$gg.causal, mode="to"),
      get.degree.in.gcc(gene.list, ggc.all.17.default.fex.ue$gg.causal, mode="from"),
      get.degree.in.gcc(gene.list, ggc.all.17.default.fex.ue$gg.causal, mode="to"),
      get.degree.in.gcc(gene.list, ggc.all.17.default.fex.ui$gg.causal, mode="from"),
      get.degree.in.gcc(gene.list, ggc.all.17.default.fex.ui$gg.causal, mode="to"))

colnames(temp.degree.ggc) <- c("ge.out", "ge.in", "gi.out", "gi.in", "ue.out", "ue.in", "ui.out", "ui.in",
                               "fge.out", "fge.in", "fgi.out", "fgi.in", "fue.out", "fue.in", "fui.out", "fui.in")

compare.between.ggc.evec <- cbind(compare.between.ggc.evec, temp.degree.ggc)

##final result##
compare.between.ggc.evec
pair.136.compare
temp.auc.array
################

#try to compare the topo sort of ggc.all.17.default.ge #
temp.graph <- graph_from_data_frame(ggc.all.17.default.ge$gg.causal, directed = T, vertices = ggc.all.17.default.ge$gene.table)
topo_sort(temp.graph, mode = "out")

#plot the bar plot 
test.histdata <- compare.between.ggc.evec[compare.between.ggc.evec$snp.freq>=7, c(1,2,3,9,18,19)]
colnames(test.histdata) <- c("gene.symbol", "out degree in EVEX", "in degree in EVEX", "snp.freq", "out degree in IGCI result", "in degree in IGCI result")
test.histdata.melt <- melt(test.histdata, id = c("gene.symbol", "snp.freq"))
colnames(test.histdata.melt)[3:4] <- c("degree", "count")
ggplot(test.histdata.melt, aes(x = gene.symbol, y=count))+
  geom_bar(aes(fill = degree), position = "dodge", stat="identity") +
  theme_classic() + xlab("Gene") + ylab("degree") 
