##infile  <- '/cluster/work/pausch/naveen/DEV/consensus/result.txt'
##outfile  <- 'stats.txt'

infile  <- snakemake@input[["infile"]]
outfile  <- snakemake@output[["outfile"]]

inf  <- read.table (infile)
apply (inf [, -c(1)],1,sum)
## ##Consersum all across all samples
mymat<- matrix (as.numeric (apply (inf [,-c(1)],2,sum)), ncol=4, byrow=T)
100*sum(diag(mymat))/sum (mymat)



concordance  <- function (mymat) {
    res  <- 100*sum (diag (mymat))/ sum (mymat)
}

mymat
NRS  <-  function (mymat) {
    nmat  <- mymat [, 2:3] ; nmat
    res  <- 100*sum ( sum(nmat [2:3,])  )  /sum (nmat) ;res
}

NRD  <- function (mymat) {
    nmat  <- mymat [1:3, 1:3]
    nmat [1,1]  <- 0
    100*    (sum(nmat) - sum (diag (nmat)) ) /sum (nmat)
}

precision  <- function (mymat) {
    nmat  <- mymat [2:3, 1:3] ; nmat
    res  <- 100 *  ( nmat [1,2] + nmat [2,3] ) /  sum (nmat)

}

k  <- 1
out  <- data.frame (matrix (nrow=nrow (inf), ncol=5)  )
colnames (out)  <- c('sample', 'conc', 'NRS', 'NRD', 'PRE')
for (k in 1:nrow (inf)) {
    out [k,1]  <- inf [k,1]
    mymat  <- as.numeric(inf [k, -c(1)])
    mymat  <- matrix (mymat,ncol=4, byrow=T)

    out [k,2]  <- concordance (mymat) 
    out [k,3]  <- NRS (mymat) 
    out [k,4]  <- NRD (mymat)
    out [k,5]  <- precision (mymat)
    
    
}
write.table (out,outfile, col.names=T,row.names=F,quote=F,sep="\t")
