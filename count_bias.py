from collections import defaultdict

##infile="/cluster/work/pausch/temp_scratch/low_pass_hap_panels/var_cal_comparison/chip_data/CONC/raw/GATK/pos_1.txt"
##outfile='counts.txt'
infile=snakemake.input.infile
outfile=snakemake.output.outfile

counts=defaultdict(dict)
snpcounts =dict ()
with open (infile, "rt") as inf:
    for lnum,line in enumerate(inf):
        spl=line.rstrip ().split ()
        if lnum==0:
            dtypes=spl[10:16]
        else:
            alleles="".join(spl[1:5])
            if len (alleles) == 4: ##keep only Bi allelic snps
                if alleles[0]==alleles[2] and alleles[1]==alleles[3]: #if the reference and alt. alleles match in the two files
                    snptype=alleles [0:2]
                    snpcounts [snptype] = snpcounts.get (snptype,0) + int(spl [7])
                    for dtype, count in zip (dtypes, spl [10:16]):
                        counts [ snptype ] [dtype] = counts [ snptype ].get (dtype,0)+ int (count)
                        
                else:
                    print (f"ref, alt alleles do not match at pos {spl[0] }")

header="\t".join (["SNP", "n",   "type", "counts"])

with open (outfile, "w") as out:
    out.write (f"{header}\n")
    for var in counts:
        for dtype in counts [var] :
            out.write (f"{var}\t{snpcounts.get(var, 0)}\t{dtype}\t{counts[var][dtype]}\n")
