# Merge bams files
# Some genotypes were sequenced multiple times. The individual runs were seperately aligned to the genome. Next, we merge the bam files together so there is one bam per genotype.
rule merge_bams:
    input:
        A = "data/interm/mark_dups/{merge_A}.dedup.bam",
        B = "data/interm/mark_dups/{merge_B}.dedup.bam"
    output:
        "data/interm/mark_dups/{geno}.{merge_A}.{merge_B}.merged.dedup.bam"
    shell:
        """
        samtools merge {output} {input.A} {input.B}
        """

# Add read groups to merged bams
rule add_rg:
    input:
        "data/interm/mark_dups/{geno}.{merge_A}.{merge_B}.merged.dedup.bam"
    output:
        bam = touch("data/interm/mark_dups/{geno}.{merge_A}.{merge_B}.merged.rg.dedup.bam")
    params:
        tmp = "/scratch/aphillip/addrg/{geno}",
        sample = "{geno}"
    run:
        shell("mkdir -p {params.tmp}")
        shell("gatk --java-options ""-Xmx4G"" AddOrReplaceReadGroups \
        -I {input} \
        -O {output.bam} \
        -RGID 4 \
        -RGLB lib1 \
        -RGPL illumina \
        -RGPU unit1 \
        -RGSM {params.sample} \
        --TMP_DIR {params.tmp} \
        --CREATE_INDEX true")
        shell("rm -rf {params.tmp}")

# Subsample high coverage bams
## subsample high coverage bams (avg cov = 29X) to 2X coverage. 2/29 = 0.06
#rule subsample:
#    input:
#        "data/interm/mark_dups/IN{high}.dedup.bam"
#    output:
#        "data/final_bams/lowcov/IN{high}.subsample.dedup.bam"
#    shell:
#        """
#        samtools view -b -s 0.06 {input} > {output}
#        """

# Quality metrics with qualimap of subsampled bams
# nr is normally 100000 and -nt is normally 8, java mem size = 48
# nw is normally 400
# for higher cov, make nr 1000 and -nt 12, java mem size = 64
# tested with nr = 10000 and nw = 400, failed
rule bamqc:
    input:
#        "data/interm/mark_dups/{geno}.{merge_A}.{merge_B}.merged.rg.dedup.bam"
        "data/final_bams/lowcov/6x_subsample/{low_geno}_{low_per}.subsample.bam",
    output:
#        "reports/bamqc/merged/{geno}.{merge_A}.{merge_B}_stats/genome_results.txt"
        "reports/bamqc/subsampled/{low_geno}_{low_per}_stats/genome_results.txt"
    params:
        dir = "reports/bamqc/subsampled/{low_geno}_{low_per}_stats"
    run: 
        shell("qualimap bamqc \
        -bam {input} \
        -nt 8 \
        -nr 1000 \
        -nw 400 \
        -outdir {params.dir} \
        -outformat HTML \
        --skip-duplicated \
        --java-mem-size=20G")    

# Subsample low coverage 6x bams to even coverage
# Goal is 0.5-1.5X coverage, only subsampling the genotypes not within that range
rule subsample_6x:
    input:
        bam_path = "data/final_bams/{low_geno}.bam"
    output:
        "data/final_bams/6x_subsample/{low_geno}_{low_per}.subsample.bam"
    params:
        frac = "{low_per}"
    shell:
        """
        samtools view -b -s {params.frac} {input.bam_path} > {output}
        """
