#' @title transfor the assay files and generate file.ped
#'
#' @description
#' \code{from.assay.to.ped} transfor the assay files and generate file.ped.
#'
#' @details
#' This function use the output of the function \code{read.assay.file}.
#' if the snp.info and sample.info file are avaliable, you should input
#' them, too. The output of this function is a ped file. The first 6 col
#' of output are about sample information. And the other col are the snp
#' from the assay files.
#'
#' @param assay output of the function \code{read.assay.file}.
#' @param snp.info a data frame containing the snp information.
#' @param sample.info a data frame containing the sample information.
#' @param is.sample.info logical. If \code{TRUE}, merge the sample information to the ped file.
#' @param col.names.convert logical. if \code{TRUE}, convert the colnames of assay to snp.ID, such as: rs001.
#' @param fam.id logical. if \code{TRUE}, the ped data will contain the family data of the samples.
#' @param is.num logical. if \code{TRUE}, the loci data in ped data will be like "1 1", "1 0" or "0 0". if \code{FALSE}, the loci data in ped data will be like "A A","A B" or "B B".
#' @param write.file logical. if \code{TRUE}, the function will write ped data down.
#' @export
#' @return a data frame with the col is sample and the row is snp.loci:
#' \item{ped}{function output}
#' @author Liangjie Liu <liuliangjie@@sjtu.edu.cn>
#' @examples
#' data.ped <- from.assay.to.ped(data.total.1536, snp.info = snp.info, sample.info = sample.info, is.sample.info = T)
#' results <- oneway(hlef ~ region, life)
#' summary(results)
#' plot(results, col="lightblue", main="Multiple Comparisons",
#'      xlab="US Region", ylab="Healthy Life Expectancy at Age 65")

###!/usr/bin/env Rscipt

###setwd("~/Desktop/data_from_HG/my_own_package/GGIGCI/")

convert.sds.to.allele <- function(rs, allele1 = NA, allele2 = NA, is.num = F) {
  if (is.num) {
    vic <- "1 1"
    fam <- "2 2"
    both <- "1 2"
  } else {
    if (is.na(allele1) | is.na(allele2)){
      stop("allele1 and allele2 are required when is.num is FALSE.\n")
    }
    vic <- paste0(allele1, " ", allele1)
    fam <- paste0(allele2, " ", allele2)
    both <- paste0(allele1, " ", allele2)
  }
  miss <- "0 0"
  out <- rep(NA, length(rs))
  for(i in 1:length(rs)) {
    n <- rs[i]
    #print(n)
    if (is.na(n)){
      out[i] <- miss
    }
    else if (n == "vic") {
      out[i] <- vic
    }
    else if (n == "fam") {
      out[i] <- fam
    }
    else if (n=="both") {
      out[i] <- both
    }
    else {
      out[i] <- miss
    }
  }
  return(out)
}

from.assay.to.ped <- function(assay, snp.info, sample.info=NA, is.sample.info = F,
                              col.names.convert = T, fam.id = F, is.num = T, write.file = F){

  if (!is.data.frame(assay) | !is.data.frame(snp.info)){
    stop("assay or snp.info must be a data.frame!\n")
  }
  if (col.names.convert) {
    if (!"ABI.ID"%in%colnames(snp.info)){
      stop("The information of ABI.ID is required when col.names.convert is TRUE.\n")
    }
  }
  if (is.sample.info) {
    if (!is.data.frame(sample.info)){
      stop("The sample.info is required when is.sample.info is TRUE.\n")
    }
  }
  if (fam.id) {
    if (!is.data.frame(sample.info)){
      stop("The sample.info is required when fam.id is TRUE.\n")
    }
  }

  if(col.names.convert) {
    snp.info[which(snp.info$ABI.ID=="" | is.null(snp.info$ABI.ID)), ]$ABI.ID <- NA
    snp.info[which(snp.info$Assays.name=="" | is.null(snp.info$Assays.name)), ]$Assays.name <- NA
    if (dim(snp.info[which(is.na(snp.info$ABI.ID) & is.na(snp.info$Assays.name)), ])[1]!=0){
      stop("Required information of ABI.ID or Assays.name missing.\n")
    }
    for(i in 1:dim(snp.info)[1]) {
      if (is.na(snp.info[i, ]$Assays.name)) {
        snp.info[i, ]$Assays.name <- snp.info[i, ]$ABI.ID
      }
      if (is.na(snp.info[i, ]$ABI.ID)) {
        snp.info[i, ]$ABI.ID <- snp.info[i, ]$Assays.name
      }
    }
    old.cnames <- snp.info$ABI.ID
    new.cnames <- snp.info$Assays.name
    temp.cname <- new.cnames

    names(temp.cname) <- old.cnames
    #print(temp.cname)
    for(i in 1:dim(assay)[2]) {
      n <- colnames(assay)[i]
      if (n%in%names(temp.cname)){
        #print(temp.cname[n])
        colnames(assay)[i] <- temp.cname[n]
      }
    }
    #print(colnames(assay))
  }

  #out <- list()
  #ped <- NA
  c <- dim(assay)[1]
  r <- dim(assay)[2]
  ped <- data.frame(Fam.id = assay$Sample.Name, Sample.Name = assay$Sample.Name,
                    Paternal.ID = rep(0, c), Maternal.ID = rep(0, c), Sex = rep(0, c),
                    Phen = rep(-9, c))
  for (i in 3:r) {
    temp.snp <- colnames(assay)[i]
    #print(temp.snp)
    #print(snp.info[which(snp.info$Assays.name==temp.snp), ]$allele1)
    allele1 <- snp.info[which(snp.info$Assays.name==temp.snp), ]$allele1
    allele2 <- snp.info[which(snp.info$Assays.name==temp.snp), ]$allele2
    #print(allele2)
    #print(length(assay[, i]))
    temp.snp.vector <- convert.sds.to.allele(tolower(assay[, i]), allele1 = allele1, allele2 = allele2, is.num = is.num)
    #print(temp.snp.vector)
    ped[, i+4] <- temp.snp.vector
    colnames(ped)[i+4] <- temp.snp
  }

  if (is.sample.info) {
    ped <- ped[, -c(1,3,4,5,6)]
    ped <- merge(sample.info, ped, by = c("Sample.Name"), all.x = T)
    ped <- ped[, c(2,1, 3:dim(ped)[2])]
  }

  return(ped)
}
