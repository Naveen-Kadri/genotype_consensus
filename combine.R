infiles  <- snakemake@input [["infiles"]]
outfile  <- snakemake@output [['outfile']]
#infiles  <- paste0('/cluster/work/pausch/naveen/DEV/consensus/CHR', 1:29, '/concordance.txt')
k  <- 1
for (k in 1:length (infiles) ) {
    inf  <- read.table (infiles [k])
    if (k ==1) {
        id = inf [,1]
        counts  <- inf [,-c(1)]
    } else {
        if (        sum (id == inf [,1]) == nrow (inf ) ) {
            counts  <- counts +inf [,-c(1)]
        }else {
            stop ('id order is not ok')
        }

    }
    cat (k, '\n')
}

out  <- data.frame (id, counts)
apply (out [,-c(1)],1,sum)
write.table (out, outfile, col.names=F,row.names=F,quote=F,sep='\t')
