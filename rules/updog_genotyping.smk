######### Filter SNPs #########
# Run diagnostics


# Hard filtering
# follow GATK best practices for hard filtering 
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2
bcftools view -i '%QUAL>=30'
