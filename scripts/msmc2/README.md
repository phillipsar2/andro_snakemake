# MSMC2 workflow
Initial scripts from JRI, modified by ARP

Notes:
- Resolution in the last 5-10Kyrs will be bad

Runs:
A v A
B v B
C v C
A v B v C

1. Anchorwave alignments generated by JGI. All-v-all alignment between haplotypes and subgenomes.
Located at: ```/group/jrigrp10/ager/data```

2. (skip - subset gvcfs instead) Make a db from the vcfs you want to combine. Use the below (changing up names as appropriate). 
This is **genomic_import.sh**

```#SBATCH --partition=high2
#SBATCH --job-name=dbimport
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=40gb
#SBATCH --time=150:00:00
#SBATCH --output=serial_test_%A_%a.log
#SBATCH --array 1-10

module load GATK/4.2.3.0

cd /group/jrigrp11/nam_gvcfs

gatk --java-options "-Xmx40g -Xms8g" GenomicsDBImport \
-V 1B.gvcf \
-V 2B.gvcf \
-V 1C.gvcf \
-V 2C.gvcf \
-V 2A.gvcf \
--genomicsdb-workspace-path ager_db_"$SLURM_ARRAY_TASK_ID"  \
--genomicsdb-vcf-buffer-size 131072000 \
--overwrite-existing-genomicsdb-workspace TRUE \
-L $SLURM_ARRAY_TASK_ID
```  

3. (skip) Then with the db made, make your new gvcf (let's call it ager.gvcf). This is **select_variants.sh**  

```#!/bin/bash
#SBATCH --partition=high2
#SBATCH --job-name=make_join_gvcf
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --time=150:00:00
#SBATCH --output=gvcf_%A_%a.log
#SBATCH --array 7

module load GATK/4.2.3.0
module load samtools

cd /group/jrigrp11/nam_gvcfs

gatk SelectVariants \
    -R Zm-B73-REFERENCE-NAM-5.0.fa \
    -V gendb://nam_db_"$SLURM_ARRAY_TASK_ID" \
    -L $SLURM_ARRAY_TASK_ID:1-70000000 \
    -L $SLURM_ARRAY_TASK_ID:90000000-180000000 \
    -O nam_"$SLURM_ARRAY_TASK_ID".vcf
```  

4. Generate file of sequenceable regions using:

zgrep ";END=" Agerardii_v1_aw.AvA.gvcf.gz| cut -f 1,2,8 | sed -e 's/ASM.*\+\;END=//g' > Agerardii_v1_aw.AvA.mask.bed

5. Generate a "clean" vcf with
(You will want to change the depth filter to whatever you want.)

perl nam_clean.pl <(zcat Agerardii_v1_aw.AvA.gvcf.gz) | grep "DP=27\|DP=25\|DP=26" | grep -v ":\.:" > Agerardii_v1_aw.AvA.clean.vcf 

**nam_clean.pl** is:

```use warnings;

my $nam=$ARGV[0];

open(NAM, '<', $nam) or die "no vcf by that name, dumbass\n";
my $skip=0;
my $line=0;
while(<NAM>){
	$line++;
	if($_=~m/^\#/){ print $_; next;}

	my @cols=split("\t",$_);

	#ballelic
	my @alt=split(",",$cols[4]);
	next unless $#alt == 1;
	my $alt_allele=""; foreach my $i (@alt){ next unless $i=~m/^[ACTG]/; $alt_allele=$i; }
	#print "$line\t$alt_allele\n";
	next if length($alt_allele)>1;

	#skip lines where ref allele >1 (deletion in some other line) and subsequent lines b/c vcf keeps including those weirdly
	my $ref=$cols[3];
	next if length($ref)>1;

	print $_;
}
```
   
Combine mask file and vcf with `zcat ager_clean.vcf.gz ager_mask.bed.gz > ager_all.txt` then sort `grep -v "#" ager_all.txt | sort  -V -k1 -k2 --parallel=10  > ager_all_sorted.txt`

Make msmc2 input file (format [here](https://github.com/stschiff/msmc-tools/blob/master/msmc-tutorial/guide.md)) `perl nam_multihep.pl ager_all_sorted.txt` 

**nam_multihep.pl** is here:
```use strict;
use warnings;

my $file=$ARGV[0];

open(FILE, '<', $file) or die "no file dufus\n";

my $space=0;
my $old_chrom="blah";
while(<FILE>){
	chomp;
	next if $_=~m/\#/;
	my @line=split("\t",$_);
	my $chrom=$line[0];
	my $start=$line[1];
	my $next=$line[2];

	if($chrom ne $old_chrom){ $space=0}

	if($next=~m/\./){
		my $ref=$line[3];
		my $alt=$line[4];
		$alt=(split(",",$alt))[0];
		print $chrom, "\t", $start, "\t", $space+1,"\t", $ref, $alt, "\n";
		$space=0;
	}
	else{
		$space=$space+$next+1-$start;
	}
	$old_chrom=$chrom;
}
```
Install msmc2. Once you have your file, run msmc2 with `~/src/msmc2/build/release/msmc2 -t 80  filename.multihetsep.txt -o filename`. You can play with the time windows (`-p` flag) if you want better resolution in some time frame, but the defaults are fine to start with.
