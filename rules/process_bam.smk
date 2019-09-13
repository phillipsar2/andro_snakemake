# Proccessing merged bam files

rule add_rg:
    input:
        config.addrg_in,
    output:
        bam = temp(touch("data/interm/addrg/{sample}.rg.bam")),
        index = temp(touch("data/interm/addrg/{sample}.rg.bai"))
    params:
        tmp = "/scratch/aphillip/addrg/{sample}"
        #sample = "{sample}",
    run:
        shell("mkdir -p {params.tmp}")
        shell("gatk AddOrReplaceReadGroups \
        -I={input} \
        -O={output.bam} \
        -RGID=4 \
        -RGLB=lib1 \
    	-RGPL=illumina \
    	-RGPU=unit1 \
    	-RGSM={{sample}} \
    	--TMP_DIR {params.tmp} \
        --CREATE_INDEX=true")
        shell("rm -rf {params.tmp}")


# This works!
rule mark_dups:
    input:
        config.mark_in
    output:
        # Make bam and index temporary
        bam = temp(touch("data/interm/mark_dups/{sample}.dedup.bam")),
        index = temp(touch("data/interm/mark_dups/{sample}.dedup.bai")),
        metrics = "qc/mark_dup/{sample}_metrics.txt"
    params:
        tmp = "/scratch/aphillip/mark_dups/{sample}"
    run:
        # Create a scratch directory
        shell("mkdir -p {params.tmp}")
        # Input bam file to output marked records. Assume bam file has been sorted. Direct to a temporary storage file (scratch).
        ## It really doesn't like my scratch files
        shell("gatk MarkDuplicates \
        -I={input} \
        -O={output.bam} \
        --METRICS_FILE={output.metrics} \
        --CREATE_INDEX=true \
        -MAX_FILE_HANDLES=1000 \
        --ASSUME_SORT_ORDER=coordinate \
        --TMP_DIR={params.tmp}")
        # Remove scratch directory
        shell("rm -rf {params.tmp}")


#gatk MarkDuplicates -I=data/sorted_reads/INCQ_PCRfree_1-McKain_252-WiDiv_43_CCGCGGTT_Andropogon_gerardii_I837_L1.sorted.bam \
#WiDiv_43_CCGCGGTT_Andropogon_gerardii_I837_L1.dedup.bam  \
#WiDiv_43_CCGCGGTT_Andropogon_gerardii_I837_L1_metrics.txt \
#--CREATE_INDEX=true \
#-MAX_FILE_HANDLES=1000 \
#--ASSUME_SORT_ORDER=coordinate \
#--TMP_DIR=/scratch/aphillip/INCQ_PCRfree_1-McKain_252-WiDiv_43_CCGCGGTT_Andropogon_gerardii_I837_L1
