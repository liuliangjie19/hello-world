#!/usr/bin/env Rscript

#setwd("~/Desktop/data_from_HG/my_own_package/GGIGCI/")

read.ped.to.snpMatrix <- function(ped, info, snp.list = NA, ped.header = T, ped.sep = "\t") {
  #selfhead(ped, nrow = 7)
  res = list()
  #require(snpStats)
  if (class(ped)!="character") {
    stop("Please, enter a valid path file for genotype data.")
  }
  file.type <- substr(ped, nchar(ped)-3, nchar(ped))
  if (file.type!=".ped") {
    stop(paste0("The file type of ", ped, " is not .ped.\n"))
  }
  if (class(info)!="character") {
    if (dim(info)[2]!=4) {
      stop("The cols of info should be Chr, Gene, Snp and Position.\n")
    }
    #res[["info"]] <- info
  } else {
    info <- read.table(info, header = T, sep = "\t")
    if (dim(info)[2]!=4){
      stop("The cols of info should be Chr, Gene, Snp and Position.\n")
    }
    #res[["info"]] <- info
  }
  #if (dim(info)[1]!=dim(ped)[1]-6) {
  #  stop("The number of SNPs in info isn't equal to the number of SNPs in the ped data.\n")
  #}
  if (ped.header) {
    temp.ped <- read.table(ped, header = ped.header, sep = ped.sep)
    if (any(is.na(snp.list))){
      snp.list <- colnames(temp.ped)[-c(1:6)]
    } else if(any(snp.list!=colnames(temp.ped[-c(1:6)]))) {
      warning("Be careful, the SNP names don't match between ped file's header and snp.list you provided.")
    }
    prefix <- format(Sys.time(), "%H%M%S%Y")
    write.table(temp.ped, file = paste0("temp", prefix, ".ped"), col.names = F, row.names = F, sep = "\t", quote = F)
    res[["ped"]] <- snpStats::read.pedfile(file = paste0("temp", prefix, ".ped"), snps = snp.list)
    #system("rm temp.ped")
    file.remove(paste0("temp", prefix, ".ped"))
  }
  else if (any(is.na(snp.list))){
    snp.list <- info[, 3]
    res[["ped"]] <- snpStats::read.pedfile(ped, snps = snp.list)
    #stop("The snp list is required when ped file without header is provided.\n")
  }
  else {
    res[["ped"]] <- snpStats::read.pedfile(ped, snps = snp.list)
  }

  #The cols of info should be Chr, Gene, Snp and Position.
  rownames(info) <- info[, 3]
  info <- info[colnames(res[["ped"]][["genotypes"]]), ]
  res[["info"]] <- info

  #res <- as(res, "GGIGCI")
  return(res)
}
