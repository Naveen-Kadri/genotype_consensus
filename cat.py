from collections import defaultdict
snpcounts=defaultdict(dict)
counts=defaultdict(dict)

#infiles = [f"/cluster/work/pausch/naveen/AUDALD/genotype_bias/COUNTS/chr_{mychr}.counts" for mychr in range (1,30) ]
#outfile='genomewide_bias.txt'


infiles=snakemake.input.infiles
outfile=snakemake.output.outfile
out=open(outfile, "w")

for mychr,myfile in enumerate(infiles):
    #print (mychr)
    with open (myfile, "rt") as inf:
        for lnum,line in enumerate(inf):
            if lnum==0:
                continue
            snp, n, type, mycount=line.rstrip().split()
            snpcounts[snp][mychr]=int (n)
            counts [snp] [type]  =   counts[snp].get (type,0)  + int (mycount)
            
            
dtypes=["0:1", "0:2", "1:0", "1:2", "2:0", "2:1"]
header = "\t".join (['SNP', 'n'] + dtypes)
out.write (f"{header}\n")
for var in counts:
    n=str (sum (list(snpcounts[var].values())))
    mydtypes =[]
    for dtype in dtypes:
        mydtypes.append( str (counts [var] [dtype]) )
    tw="\t".join ([var, n] + mydtypes)
    out.write (f"{tw}\n")

out.close ()

