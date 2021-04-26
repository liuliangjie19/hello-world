require(ggplot2)
require(openssl)
require(gdata)
getwd()
source("R/read.assay.file.R")
data.total.1536 <- read.assay.file(dir = "~/Desktop/data_from_HG/data_first_20190116/", pattern = "AZ")
#selfhead(data.total.1536)
#data.total.1536
ls()
head(snp2gene)
dim(snp2gene)
??read.xls
require(XRJulia)
require(JuliaCall)
fun1 <- juliaEval("
                  function f(x)
                   sum = 0
                   for xi in x
                    sum = sum + xi
                   end
                   return(sum)
                  end")
fun1.julia <- JuliaFunction(fun1)
fun1.julia(snp.info.gene1)
fun1.julia(phen)
length(phen)
rm(fun1.julia)
rm(fun1)
funtest <- juliaEval("
                     function statcalculate(phen)
                       b = [1,2]
                       b[3] = 3
                       return(b)
                     end")
funtest <- JuliaFunction(funtest)
funtest(xt1)
