from collections import defaultdict
import gzip

truth_file = "/cluster/work/pausch/naveen/SNPDATA/hd/CHR25/phased.vcf.gz"
test_file=  "/cluster/work/pausch/naveen/VARCAL_2022_03/CHR25/PHASED/cohort_25_beagle4.1.vcf.gz"
out_file = 'indwise.txt'

#truth_file=snakemake.input.truth_file
#test_file = snakemake.input.test_file
#out_file = snakemake.output.out_file

truth = gzip.open (truth_file,  'rt')
test = gzip.open (test_file, 'rt')
out= open (out_file, 'w')

print ("checking for overlap")
for line in truth:
    if line [0:6] == "#CHROM":
        true_ids=line.strip().split() [9:]
        break
for line in test:
    if line [0:6] == "#CHROM":
        test_ids=line.strip().split() [9:]
        break
overlapping_ids=set(true_ids).intersection (set(test_ids))
print (f"Number of samples in the first & second file {len (true_ids):,} &  {len (test_ids):,}")
print (f"Number of overlapping ids : {len (overlapping_ids) } " )

if len(overlapping_ids) ==0 :
    print ("No overlapping id in the two files ; exiting the  program")
    exit ()


recode={
    "0/0" : "0",    "0/1" : "1",    "1/0" : "1",    "1/1" : "2", "./.":"5",
    "0|0" : "0",    "0|1" : "1",    "1|0" : "1",    "1|1" : "2", ".|.":"5"
}
##order to write genotype comparison [to read as matrix]
recoded = ['0', '1', '2', '5']
write_order = []
for i in recoded:
    for j in recoded:
        write_order.append (i + ":" + j)

result = defaultdict (dict)
nov = 0
nc = 0
def GTcompare (true_ids, true_spl, test_ids, test_spl):
    true_gts = [recode.get (el[0:3], 'NA') for el in true_spl [9:] ]
    test_gts = [recode.get (el[0:3], 'NA') for el in test_spl[9:] ]
    x=sum (list (map (int, true_gts)))
    print (true_gts)
    print (x)
    exit ()

    true_geno = dict (zip (true_ids, true_gts))
    test_geno = dict (zip (test_ids, test_gts))
    for my_testid, my_testgt in test_geno.items ():
        if my_testid in true_geno:
            gcomb = ":".join ( [test_geno [my_testid], true_geno[my_testid] ] )
            result [my_testid] [gcomb] = result [my_testid].get (gcomb,0) +1
            
test_pos=None
ovpos = open ('ovpos.txt' ,'w')
for line in truth:
    true_spl=line.rstrip ().split ()
    true_pos = int (true_spl [1])
    #if true_pos > 5_000_000:
        #break
    if test_pos is not None:
        if true_pos ==test_pos:
            nov +=1
            ovpos.write (f'{nov}\tb1\t{test_pos}\n')
            print (f'nov = {nov} at position {test_pos:,}')
            GTcompare (true_ids, true_spl, test_ids, test_spl)
        elif true_pos < test_pos:
            continue
    for line in test:
        test_spl=line.rstrip().split()
        test_pos = int(test_spl[1])
        if test_pos == true_pos:
            ovpos.write (f'{nov}\tb2\t{test_pos}\n')
            nov +=1
            print (f'nov = {nov} at position {test_pos:,}')
            GTcompare (true_ids, true_spl, test_ids, test_spl)
        elif test_pos > true_pos:
            break
    
ovpos.close()
truth.close ()
test.close ()

allresult =dict ()
for myid in result:
    out.write (f'{myid}')
    for gcomb in write_order:
        out.write (f'\t{result[myid].get (gcomb,0) }')
        allresult [gcomb] = allresult.get (gcomb,0) + result[myid].get (gcomb,0)
    out.write ('\n')
out.close()
print (f'Number of overlapping positions : {nov}')


