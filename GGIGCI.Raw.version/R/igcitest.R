for(file in dir(path, pattern = "AZ")){
  #print(file)
  temp.data <- read.table(paste0(path, file), skip = 11, sep = "\t", header = T)
  temp.snp.name <- strsplit(file, split = "-")[[1]][3]
  temp.plate <- strsplit(file, split = "-")[[1]][4]
  #temp.Well <- temp.data$Well
  #temp.Call <- temp.data$Call
  #temp.Name <- temp.data$Sample.Name
  temp.date.frame <- temp.data[, c(1,2,6)]
  
}

for(file in dir(path , pattern = "plate2")) {
  #print(file)
  temp.data <- read.table(paste0(path, file), skip = 11, sep = "\t", header = T)
  #print(dim(temp.data))
  if (dim(temp.data)[1]!=384) {
    print(file)
    print(dim(temp.data))
    print(head(temp.data))
  }
}

file <- "./data_first_20190116/AZ-SLC6A4-C_7911132_10-plate2-20051010-001982.txt"
temp.data <- read.table(file, skip = 11, header = T, sep = "\t")
temp.data[which(row.names(temp.data)!=temp.data$Well),]

for(i in 1:dim(snp.info)[1]) {
  x <- snp.info$ABI.ID[i]
  
  snp.info$ABI.ID[i] <- gsub("_+", "_", x)
  print(x)
}
snp.info$ABI.ID


geno_data <- read.table("data_first_20190116/data.total", header = T, sep = "\t")
snp2gene <- read.table("data_first_20190116/data.txt", header = T, sep = "\t")
snp.map <- read.table("data_first_20190116/data.info", header = F, sep = "\t")

gene.list <- unique(snp2gene$Genenames)
#strsplit(snp.map$V2, split = ":")
snp.map <- data.frame(chr = snp2gene$Chromosome, snp_id = snp2gene$SNPnames, distant = rep(0, dim(snp2gene)[1]), posi = snp2gene$Position)

names(table(geno_data[,10]))
temp.geno_data <- geno_data

snp.pos <- snp.map$posi
names(snp.pos) <- snp.map$snp_id

for(i in 7:98) {
  t.name <- names(table(geno_data[, i]))
  if (length(t.name)!=4){
    message(colnames(geno_data)[i], "wrong!")
    print(t.name)
  }
  geno_data[which(geno_data[, i]==t.name[1]), i] <- -9
  geno_data[which(geno_data[, i]==t.name[2]), i] <- 0
  geno_data[which(geno_data[, i]==t.name[3]), i] <- 1
  geno_data[which(geno_data[, i]==t.name[4]), i] <- 2
}
#temp.geno_data[1, ]

geno.num.gene.stat <- geno_data[, c(1:6)]
geno.fex.data <- geno.num.gene.stat
gene.n.basis <- rep(0, length(gene.list))
names(gene.n.basis) <- gene.list

for(gene in gene.list) {
  print(gene)
  snp.in.gene = snp2gene[snp2gene$Genenames%in%gene,]$SNPnames
  #print(gene)
  geno.num.gene = cbind(geno_data[,c(1:6)], geno_data[,colnames(geno_data)%in%snp.in.gene])
  geno.num.gene[geno.num.gene<0]=NA
  geno.num.gene = na.omit(geno.num.gene)
  geno.num.gene.stat = geno.num.gene[,c(1:6)]
  geno.num.gene = geno.num.gene[,-c(1:6)]
  print(dim(geno.num.gene))
  nsnps = dim(geno.num.gene)[2]
  gil = list()
  pos = as.numeric(snp.pos[names(snp.pos)%in%colnames(geno.num.gene)])
  if(!is.null(nsnps)) {
    if ( is.null(pos) ) {
      pos <-(0:(nsnps-1))/(nsnps-1)
    }
    else {
      idx = order(pos)
      geno.num.gene = geno.num.gene[,idx]
      pos = pos[idx]
      pos<-(pos-pos[1])/(pos[nsnps]-pos[1])
      print(pos)
    }
    # use PCA to determine the number of basis
    geno.num.gene <- apply(geno.num.gene, 2, as.numeric)
    eigenval<-prcomp(geno.num.gene)$sd^2
    #print(gene)
    sum_eigen=sum(eigenval)
    tmp=0
    n_of_basis=0
    for(i in 1:length(eigenval))
    {
      tmp=eigenval[i]+tmp
      n_of_basis=i;
      if(tmp>=0.8*sum_eigen) {
        break
      }
    }
    n_of_basis=floor(n_of_basis/2)*2+1 
    print(n_of_basis)
    gene.n.basis[gene] <- n_of_basis
    #make n_of_basis_A the odd number
    #end of setting the basis number
    frange <-c(pos[1], pos[length(pos)])
    fbasis<-create.fourier.basis(frange,nbasis=n_of_basis)
    phi=eval.basis(pos,fbasis);
    gil$zeta = t(ginv(t(phi)%*%phi)%*%t(phi)%*%t(geno.num.gene))
    gil$zeta = cbind(geno.num.gene.stat$Sample.Name, gil$zeta)
    colnames(gil$zeta) <- c("Sample.Name", paste0(gene, c(1:n_of_basis)))
  }
  else if (is.null(nsnps)) {
    gene.n.basis[gene] <- 1
    gil$zeta <- cbind(geno.num.gene.stat$Sample.Name, geno.num.gene)
    colnames(gil$zeta) <- c("Sample.Name", paste0(gene, "1"))
  }
  geno.fex.data <- merge(geno.fex.data, gil$zeta, by = c("Sample.Name"), all.x = T)
}

inter <- combn(gene.list, 2)
gene.graph <- data.frame(gene1 = rep(NA, ncol(inter)), gene2 = rep(NA, ncol(inter)), 
                         S = rep(0, ncol(inter)), D = rep(0, ncol(inter)))
for(i in 1:ncol(inter)){
  gene1 <- inter[1, i]
  gene2 <- inter[2, i]
  n1 <- gene.n.basis[gene1]
  n2 <- gene.n.basis[gene2]
  v1 <- paste0(gene1, c(1:n1))
  v2 <- paste0(gene2, c(1:n2))
  #if (n1==1){
  #  v1 <- c(gene1)
  #}
  #if (n2==1){
  #  v2 <- c(gene2)
  #}
  d1 <- as.data.frame(geno.fex.data[, v1])
  d2 <- as.data.frame(geno.fex.data[, v2])
  s <- 0
  for(k in 1:n1){
    for (j in 1:n2){
      temp <- cbind(d1[,k], d2[,j])
      temp <- na.omit(temp)
      #print(igci(abs(as.numeric(temp[,1])), abs(as.numeric(temp[,2]))))
      s <- s+igci(abs(as.numeric(temp[,1])), abs(as.numeric(temp[,2])))
    }
  }
  print(paste0(gene1, " to ", gene2, ":", as.character(s)))
  if (s>0) {
    gene.graph[i, ] <- c(gene2, gene1, s, 1)
  }else {
    gene.graph[i, ] <- c(gene1, gene2, s, s/abs(s))
  }
}

require(igraph)
gene.v <- data.frame(gene.list, gene.n.basis)

g <- graph_from_data_frame(gene.graph, directed = T, vertices = gene.v)
plot(g)
