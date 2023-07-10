import gzip
import sys
import numpy as np


#file1="/cluster/work/pausch/naveen/BOVREG/1KCLEAN/BREEDWISE/COMBINED/CHR25/combined.vcf.gz"
#file2="/cluster/work/pausch/naveen/ALLIMPUTE/hd/CHR25/ref_alt_beagle5.vcf.gz"
#file2="/cluster/work/pausch/naveen/BOVREG/1KCLEAN/CHR25/BEAGLE/cleaned_beagle4.1.vcf.gz"


file1="/cluster/work/pausch/naveen/ALLIMPUTE/hd/CHR25/ref_alt_beagle5.vcf.gz"
file2="/cluster/work/pausch/naveen/BOVREG/1KCLEAN/BREEDWISE/COMBINED/CHR25/combined.vcf.gz"


print (f"First file   : {file1}")
print (f"Second  file : {file2}")
f=gzip.open (file1, "rt") 
s=gzip.open (file2, "rt") 

print ("checking for overlap")
for line in f:
    if line [0:6] == "#CHROM":
        spl=line.strip().split()
        all_fids = np.array(spl [9:])
        break
for line in s:
    if line [0:6] == "#CHROM":
        spl=line.strip().split()
        all_sids = np.array(spl [9:])
        break
overlapping_ids=np.intersect1d(all_fids, all_sids)

##indices for the overlapping_ids in the two files
f_select=np.isin (all_fids, overlapping_ids)
s_select=np.isin (all_sids, overlapping_ids)
fids=all_fids [f_select]
sids=all_sids [s_select]

#find a sorting index to sort the genotypes of the second file
#then the order of the ids will be as the "fids"
sorting_index=[] ##for the second file to follow the order of fids !
for fi in fids:
    for i,si in enumerate (sids):
        if fi==si:
            sorting_index.append (i)
            break

print (f"Number of samples in the first file {len (all_fids):,}")
print (f"Number of samples in the second file {len (all_sids):,}")
print (f"Number of overlapping ids : {len (overlapping_ids) } " )
if len(overlapping_ids) ==0 :
    print ("No overlapping id in the two files ; exiting the  program")
    exit ()
    
spos=None
recode={
    "0/0" : 0,    "0/1" : 1,    "1/0" : 1,    "1/1" : 2,
    "0|0" : 0,    "0|1" : 1,    "1|0" : 1,    "1|1" : 2
}

n_pos_overlap=0
per_id=np.zeros (len (overlapping_ids))
per_id_ncomp=np.zeros (len (overlapping_ids))
out1=open("consensus_pos.txt", "w")
out1.write ("Position\tFile1Ref\tFile1Alt\tFile2Ref\tFile2Alt\tFreqF\tFreqS\tncompared\tnmatch\tfraction\n")
out2=open("consensus_id.txt", "w")
out2.write ("id\tncompared\tnmatched\tfraction\n")


test =open ("test.txt","w")
def GTcompare (): #fgts, fids, sgts and sids
        
        
    global per_id
    global per_id_ncomp
    global n_pos_overlap
    n_pos_overlap+=1
    #freq
    f_freq= np.mean(fgts[fgts >=0])/2
    s_freq= np.mean(sgts[sgts >=0])/2
    
    sums=fgts + sgts
    n = sum(sums >= 0)
    per_id_ncomp=per_id_ncomp + (sums >=0 )
    matches=fgts == sgts
    #at the missing site the genotypes will match if both are missing ! so this needs to be handled! easy thing woould be code them differently in the two files
    per_id=per_id+ matches
    nmatch=sum (matches)
    out1.write (f"{fpos}\t{fref}\t{falt}\t{sref}\t{salt}\t{f_freq}\t{s_freq}\t  {n}\t{nmatch}\t{nmatch/n}\n")
    #if n_pos_overlap % 100 == 0:
    print (f"number of overlapping loci : {n_pos_overlap:,} at {fpos:,}")
 
        
for fline in f:
    fspl=fline.rstrip ().split ()
    fpos = int (fspl [1])
    fref = fspl [3]
    falt = fspl [4]
    fgts=np.array (fspl[9:])
    ##keep genotypes only for the overlapping ids
    fgts=fgts [f_select]
    fgts=np.array ([recode.get (gt[0:3]  ,-7)  for gt in fgts ],dtype="int8")
    #here
    ##if fpos > 100000:
        ##break
    
    if spos is not None:
        if fpos ==spos:
            GTcompare ()
        if spos > fpos:
            continue  
    for sline in s:
        sspl=sline.rstrip().split()
        spos = int(sspl[1])
        sref = sspl [3]
        salt = sspl [4]
        sgts=np.array(sspl[9:])
        sgts =sgts [s_select]
        sgts = sgts [sorting_index]
        sgts = np.array ([recode.get (gt[0:3],-9) for gt in sgts],dtype="int8")

        if spos > fpos:
            break
        
        if fpos == spos:
            GTcompare ()

f.close ()
s.close ()
print (f"Number of overlapping positions : {n_pos_overlap:,}")

for myid, ncomp, matched  in zip(fids,  per_id_ncomp, per_id):
    out2.write (f"{myid}\t{ncomp}\t{matched}\t{matched/ncomp}\n")


out1.close ()
out2.close ()





