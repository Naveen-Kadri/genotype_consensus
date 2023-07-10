infiles = ["/cluster/work/pausch/temp_scratch/low_pass_hap_panels/var_cal_comparison/chip_data/CONC/trim_seq/GATK/raw/pos_1.txt", "/cluster/work/pausch/temp_scratch/low_pass_hap_panels/var_cal_comparison/chip_data/CONC/trim_seq/DV/raw/pos_1.txt"]
outfile = "todel_common.txt"

out =open (outfile, "w")
snps = dict ()
nfirst = 0
nsecond = 0
ncommon = 0
for fnum, myfile in enumerate (infiles):
    with open (myfile, "rt") as inf:
        for lnum,line in enumerate(inf):
            if lnum ==0:
                continue
            pos=line.split ("\t") [0]
            if fnum == 0:
                snps [pos] =1
                nfirst +=1
            else:
                nsecond +=1
                if pos in snps:
                    ncommon+=1
                    out.write (f"{pos}\n")
out.close ()
print (f"Number of snps in the first file, second file & common in the two  : {nfirst}, {nsecond} & {ncommon}")

