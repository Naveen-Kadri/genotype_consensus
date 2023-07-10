localrules:combine
chromosomes = range (1,30)
#truth_file = '/cluster/work/pausch/to_share/audald/conc_test/multi_sample/truth_set.vcf.gz'
#test_file = '/cluster/work/pausch/to_share/audald/conc_test/multi_sample/query_GLIMPSE.vcf.gz'

test_file = '/cluster/work/pausch/temp_scratch/low_pass_hap_panels/GLIMPSE/BSW_panel/full_cohort/Two/imputed/{chr}.vcf.gz'
truth_file = '/cluster/work/pausch/temp_scratch/low_pass_hap_panels/truth_set/var_cal/imputation/{chr}_beagle.vcf.gz'
OUT_DIR = '/cluster/work/pausch/naveen/DEV/consensus'
rule all:
    input:
        OUT_DIR + '/stats.txt'

rule compare_vcf:
    input:
        truth_file=truth_file,
        test_file=test_file
    output:
        out_file=OUT_DIR + '/CHR{chr}/concordance.txt'
    resources:
        mem_mb=4000,
        walltime= '04:00'
    script:
        "concordance_v5.py"

rule combine:
    input:
        infiles = expand (OUT_DIR + '/CHR{chr}/concordance.txt', chr=chromosomes)
    output:
        outfile =OUT_DIR  + '/result.txt'
    script:
        'combine.R'    


rule stats:
    input:
        infile=rules.combine.output.outfile
    output:
        outfile =OUT_DIR + '/stats.txt'
    script:
        'stats.R'
