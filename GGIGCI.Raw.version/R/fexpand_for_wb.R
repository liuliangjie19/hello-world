library(FRGEpistasis)
library(fda)
library(MASS)
library(stats)
library(BiocManager)
BiocManager::install("factoextra")
library(factoextra)

library(openxlsx)
snplist <- read.csv("snp_map.txt", header = F, sep = "\t")
snpinfo <- read.xlsx("HSCR数据整理-原始(58个SNPs).xlsx", 
                     startRow = 2, colNames = TRUE,sheet = 2)
snpinfo <- snpinfo[-c(59,60,61,62), ]
snpinfo$Gene[1:14] <- snpinfo$Gene[1]
snpinfo$Gene[15:38] <- snpinfo$Gene[15]
snpinfo$Gene[39:49] <- snpinfo$Gene[39]
snpinfo$Gene[50:58] <- snpinfo$Gene[50]
genelist <- unique(snpinfo$Gene)
#genelist <- genelist[1:4]
gene_chr <- c(1,3,3,10)
gene_start <- c(116242624, 27525911, 171318195, 43187207)
gene_end <- c(116311426, 27757440, 171528284, 43248359)
gene_map <- cbind(genelist, gene_chr, gene_start, gene_end)
colnames(gene_map) <- c("Gene_Symbol", "Chromosome", "Start", "End")
gene_map <- as.data.frame(gene_map)
snp_map <- cbind(snpinfo$CHR, snpinfo$SNP, rep(0, dim(snpinfo)[1]), snpinfo$POS)
colnames(snp_map) <- c("chr", "snp_id", "distant", "posi")
snp_map <-as.data.frame(snp_map)
raw.data <- read.xlsx("HSCR数据整理-原始(58个SNPs).xlsx", sheet = 1,
                      colNames = TRUE)
library(synbreed)
ge <- raw.data[,7:64]
info <- raw.data[,1:6]
ge[ge == "0 0"] <- NA
ge <- create.gpData(geno=ge)
ge.digit <- codeGeno(gpData = ge,label.heter = "alleleCoding", maf = 0.01, nmiss = 0.1,
                   impute = TRUE, impute.type = "random", verbose = TRUE)
ge.digit <- ge.digit$geno
ge.data <- cbind(info, ge.digit)
rng = 0
snp.list <- colnames(ge.data)[-c(1:6)]
snp.id <- snp_map$snp_id
snp.pos <- as.vector(snp_map$posi)
names(snp.pos) <- as.vector(snp.id)

for(gene1 in 1:length(genelist)) {
  for(gene2 in gene1:length(genelist)) {
    if (gene1!=gene2){
      print(genelist[gene1])
      print(genelist[gene2])
      snp.in.gene1 <- snpinfo[snpinfo$Gene%in%genelist[gene1],]$SNP
      snp.in.gene2 <- snpinfo[snpinfo$Gene%in%genelist[gene2],]$SNP
      gene1.snp <- ge.digit[,colnames(ge.digit)%in%snp.in.gene1]
      gene2.snp <- ge.digit[,colnames(ge.digit)%in%snp.in.gene2]
      phen <- info$`1=control,.2=case`-1
      print(Aggregator(Y = phen, gene1.snp, gene2.snp))
    }
  }
}

for(gene in genelist) {
  print(gene)
  snp.in.gene <- snpinfo[snpinfo$Gene%in%gene,]$SNP
  gene.snp <- ge.digit[,colnames(ge.digit)%in%snp.in.gene]
  print(dim(gene.snp))
  nsnps <- dim(gene.snp)[2]
  gil <- list()
  pos = as.numeric(snp.pos[names(snp.pos)%in%colnames(gene.snp)])
  if(nsnps>1) {
    if ( is.null(pos) ) {
      pos <-(0:(nsnps-1))/(nsnps-1)
    }
    else {
      idx = order(pos)
      gene.snp = gene.snp[,idx]
      pos = pos[idx]
      pos<-(pos-pos[1])/(pos[nsnps]-pos[1])
      print(pos)
    }
    # use PCA to determine the number of basis
    eigenval<-prcomp(gene.snp)$sd^2
    sum_eigen=sum(eigenval)
    tmp=0
    n_of_basis=0
    for(i in 1:length(eigenval)){
      tmp=eigenval[i]+tmp
      n_of_basis=i;
      if(tmp>=0.8*sum_eigen) {
        break
      }
    }
    n_of_basis=floor(n_of_basis/2)*2+1 
    print(n_of_basis)
    #make n_of_basis_A the odd number
    #end of setting the basis number
    frange <-c(pos[1], pos[length(pos)])
    fbasis<-create.fourier.basis(frange,nbasis=n_of_basis)
    phi=eval.basis(pos,fbasis);
    gil$zeta = t(ginv(t(phi)%*%phi)%*%t(phi)%*%t(gene.snp))
    gil$zeta = cbind(info, gil$zeta)
  }
  else {
    gil$zeta  <-cbind(info, gene.info)
  }
  filename = paste(gene, "fex.csv", sep = "_")
  write.csv(gil$zeta, file = filename, quote = F)
}
print(genelist)

pval.mat = matrix(data = NA, nrow = length(genelist), ncol = length(genelist))
rownames(pval.mat) = colnames(pval.mat) = genelist

for(gene1 in 1:length(genelist)) {
  gene.fex1 = read.csv(paste(genelist[gene1], "fex.csv", sep = "_"), header = FALSE, skip = 1)
  n1 = dim(gene.fex1)[2]-7
  #print(n1)
  for(gene2 in 1:length(genelist)) {
    if (gene2 == gene1){
      pval.mat[gene1, gene2] = 1
      next
    }
    gene.fex2 = read.csv(paste(genelist[gene2], "fex.csv", sep = "_"), header = FALSE, skip = 1)
    n2 = dim(gene.fex2)[2]-7
    #gene.fex = merge(gene.fex1, gene.fex2, by=c("X", "Sample.ID", "Year.of.birth", "Sex", "Age.of.onset","age", "Stat"), all=FALSE)
    gene.fex <- merge(gene.fex1, gene.fex2, by=c("V1", "V2", "V3", "V4", "V5", "V6", "V7"))
    gene.fex$V7 <- gene.fex$V7 - 1
    pval.vec = rep(NA, times = n1*n2)
    n = 1
    for(i in 8:(8+n1-1)) {
      for(j in (7+n1+1):(7+n1+n2)) {
        #print(anova(glm(gene.fex[,2]~gene.fex[,i]*gene.fex[,j], family="binomial"), test="Chisq"))
        pval.vec[n] = anova(glm(gene.fex[,7]~gene.fex[,i]*gene.fex[,j], family="binomial"), test="Chisq")[4,5]
        n = n+1
      }
    }
    print(gene1)
    print(gene2)
    print(min(pval.vec))
    pval.mat[gene1, gene2] = pval.mat[gene2, gene1] = min(pval.vec)
  }
}

write.csv(pval.mat, file = "pvalmatrix.csv")


pvalmat1 = read.csv("result_p_aggregator.csv", header = T)
pvalmat1 = rbind(pvalmat1, rep(NA, times = 17))
colnames(pvalmat1)[6] = "G72/G30.txt"
rownames(pvalmat1) = colnames(pvalmat1)
colnames(pvalmat1)
rownames(pvalmat1)
for(i in 1:17){
  for(j in i:17){
    if (i==j){
      pvalmat1[j,i]=1
    }
    pvalmat1[j,i] = pvalmat1[i,j]
  }
}
#pval.mat[rownames(pval.mat)%in%gene_list1[1],colnames(pval.mat)%in%gene_list1[2]]
#pvalmat1[rownames(pvalmat1)%in%paste(gene_list1[1], "txt", sep="."), colnames(pvalmat1)%in%paste(gene_list1[1], "txt", sep=".")]
write.csv(pvalmat1, file = "pvalmatrixbefore.csv")
pcomp = matrix(data = NA, nrow = 2, ncol = length(gene_list1)*length(gene_list1))
comp = 1
p1 = 0
p2 = 0
for(gene1 in 1:length(gene_list1)) {
  for(gene2 in 1:length(gene_list1)){
    pval1 = pval.mat[rownames(pval.mat)%in%gene_list1[gene1],colnames(pval.mat)%in%gene_list1[gene2]]
    pval2 = pvalmat1[rownames(pvalmat1)%in%paste(gene_list1[gene1], "txt", sep="."), colnames(pvalmat1)%in%paste(gene_list1[gene2], "txt", sep=".")]
    #print(gene_list1[gene1])
    #print(gene_list1[gene2])
    #print(pval1)
    #print(pval2)
    if (pval1<0.05){
      print(paste("1", gene_list1[gene1], gene_list1[gene2], pval1, pval2, sep = "__"))
      p1 = p1 + 1
    }
    if (pval2<0.05){
      print(paste("2", gene_list1[gene1], gene_list1[gene2], pval1, pval2, sep = "__"))
      p2 = p2 + 2 
    }
    pcomp[1,comp] = pval1
    pcomp[2,comp] = pval2
    comp = comp+1
  }
  #break
}
plot(pcomp[1,],pcomp[2,])
x = c(0:1)
y = x
lines(y~x)
pcomp1=pcomp[,which(pcomp[1,]<0.05 | pcomp[2,]<0.05)]
plot(pcomp1[1,], pcomp1[2,])
lines(y~x)
plot(pcomp2[1,], pcomp2[2,])
lines(y~x)
plot(density(COMT.fex[which(COMT.fex$Stat==1),]$X1),xlim=c(-1,3), ylim=c(0,2), col=1)
lines(density(COMT.fex[which(COMT.fex$Stat==0),]$X1), col=2)
lines(density(COMT.fex[which(COMT.fex$Stat==1),]$X2), col=3)
lines(density(COMT.fex[which(COMT.fex$Stat==0),]$X2), col=4)
lines(density(COMT.fex[which(COMT.fex$Stat==1),]$X3), col=5)
lines(density(COMT.fex[which(COMT.fex$Stat==0),]$X3), col=6)
lines(density(COMT.fex[which(COMT.fex$Stat==1),]$X4), col=7)
lines(density(COMT.fex[which(COMT.fex$Stat==0),]$X4), col=8)
lines(density(COMT.fex[which(COMT.fex$Stat==1),]$X5), col=9)
lines(density(COMT.fex[which(COMT.fex$Stat==0),]$X5), col=10)
lines(density(COMT.fex[which(COMT.fex$Stat==1),]$X4))
lines(density(COMT.fex[which(COMT.fex$Stat==1),]$X5))
