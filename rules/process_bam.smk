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
        -I={input} \
        -O={output.bam} \
        -RGID=4 \
        -RGLB=lib1 \
        -RGPL=illumina \
        -RGPU=unit1 \
        -RGSM={params.sample} \
        --TMP_DIR {params.tmp} \
        --CREATE_INDEX=true")
        shell("rm -rf {params.tmp}")
