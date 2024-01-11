farm🌿:~$ cd /group/jrigrp11/
farm🌿:jrigrp11$ ls
aphillip    contact  nam10.log  pmagalang           single_read  tga1              zeasyn
assemblies  cstark   ndhamra2   relate_data_regina  sodell       zeamap_diversity
farm🌿:jrigrp11$ cd aphillip/
farm🌿:aphillip$ ls
andropogon  cg_andropogon
farm🌿:aphillip$ cd andropogon/
farm🌿:andropogon$ ls
-                          slurm-5472050.out  slurm-5472079.out  slurm-5472154.out  slurm-5473446.out  slurm-5475944.out
6x.lowcov.samples.csv      slurm-5472051.out  slurm-5472080.out  slurm-5472155.out  slurm-5473447.out  slurm-5475945.out
6x_bams_to_subset.tsv      slurm-5472052.out  slurm-5472081.out  slurm-5472156.out  slurm-5473448.out  slurm-5475946.out
6x_lowcov_bams.txt         slurm-5472053.out  slurm-5472082.out  slurm-5472157.out  slurm-5473449.out  slurm-5475947.out
9x.lowcov.samples.csv      slurm-5472054.out  slurm-5472083.out  slurm-5472158.out  slurm-5473450.out  slurm-5475948.out
HAP1_A.gene_exons.gff3     slurm-5472055.out  slurm-5472084.out  slurm-5472159.out  slurm-5473451.out  slurm-5475949.out
README.md                  slurm-5472056.out  slurm-5472085.out  slurm-5472160.out  slurm-5473452.out  slurm-5475950.out
R_libs                     slurm-5472057.out  slurm-5472086.out  slurm-5472161.out  slurm-5473453.out  slurm-5475951.out
Rplots.pdf                 slurm-5472058.out  slurm-5472087.out  slurm-5472162.out  slurm-5473454.out  slurm-5475952.out
Snakefile                  slurm-5472059.out  slurm-5472134.out  slurm-5472163.out  slurm-5473455.out  slurm-5475953.out
__pycache__                slurm-5472060.out  slurm-5472135.out  slurm-5472164.out  slurm-5473456.out  slurm-5475954.out
allgenos_w_ploidy.csv      slurm-5472061.out  slurm-5472136.out  slurm-5472165.out  slurm-5473457.out  slurm-5475955.out
bams_nquire.csv            slurm-5472062.out  slurm-5472137.out  slurm-5472166.out  slurm-5473458.out  slurm-5475956.out
bams_to_merge.tsv          slurm-5472063.out  slurm-5472138.out  slurm-5472167.out  slurm-5473459.out  slurm-5475957.out
config.py                  slurm-5472064.out  slurm-5472139.out  slurm-5472168.out  slurm-5473460.out  slurm-5475958.out
data                       slurm-5472065.out  slurm-5472140.out  slurm-5472169.out  slurm-5473461.out  slurm-5475959.out
environment.ym             slurm-5472066.out  slurm-5472141.out  slurm-5472170.out  slurm-5473462.out  slurm-5475960.out
fastp.html                 slurm-5472067.out  slurm-5472142.out  slurm-5472171.out  slurm-5473463.out  slurm-5475961.out
highcov                    slurm-5472068.out  slurm-5472143.out  slurm-5473435.out  slurm-5473464.out  slurm-5475962.out
highcov_subsample_map.txt  slurm-5472069.out  slurm-5472144.out  slurm-5473436.out  slurm-5473465.out  slurm_log
highcov_w_ploidy.csv       slurm-5472070.out  slurm-5472145.out  slurm-5473437.out  slurm-5473466.out  smudgeplot
logs                       slurm-5472071.out  slurm-5472146.out  slurm-5473438.out  slurm-5473467.out  status.py
models                     slurm-5472072.out  slurm-5472147.out  slurm-5473439.out  slurm-5473468.out  structure
multidog.jpeg              slurm-5472073.out  slurm-5472148.out  slurm-5473440.out  slurm-5473469.out  submit.json
notebooks                  slurm-5472074.out  slurm-5472149.out  slurm-5473441.out  slurm-5473470.out  submit.sh
qc                         slurm-5472075.out  slurm-5472150.out  slurm-5473442.out  slurm-5473471.out  treemix
reports                    slurm-5472076.out  slurm-5472151.out  slurm-5473443.out  slurm-5473472.out  vcf
rules                      slurm-5472077.out  slurm-5472152.out  slurm-5473444.out  slurm-5475942.out  workshop
scripts                    slurm-5472078.out  slurm-5472153.out  slurm-5473445.out  slurm-5475943.out
farm🌿:andropogon$ cd scripts/
farm🌿:scripts$ l
l: command not found
sfarm🌿:scripts$ sls
Command 'sls' not found, but there are 19 similar ones.
farm🌿:scripts$ ls
Entropy.sh   archive              div_test.sh     genoDPfilter.sh    glmat2mpgl.R    inbreeding.sh  ngsf.sh      strict_missing_data.R
README.md    assessconvergence.R  fastqc.sh       genome_minimap.sh  gltable2mpgl.R  jellyfish.sh   ploidFIle    vcf2ADmatrix.R
angsdpca.sh  bwa_index.sh         genoDPfilter.R  get_af.R           hmmploidy.sh    mapping3.sh    plotadmix.R
farm🌿:scripts$ cd ..
farm🌿:andropogon$ cd rules/
farm🌿:rules$ ls
aneuploidy.smk  calling.smk           filtering.smk  old_rules      pop_struc_highcov.smk  scripts
angsd.smk       determine_ploidy.smk  mapping.smk    pop_struc.smk  process_bam.smk
farm🌿:rules$ cd ..
farm🌿:andropogon$ cd scripts/
farm🌿:scripts$ mkdir msmc2
farm🌿:scripts$ cd msmc2/
farm🌿:msmc2$ ls
farm🌿:msmc2$ Make a db from the vcfs you want to combine. Use the below (changing up names as appropriate). This is **genomic_import.sh**

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

Then with the db made, make your new gvcf (let's call it ager.gvcf). This is **select_variants.sh**  

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

Generate file of sequenceable regions using `grep ";END=" ager.gvcf| cut -f 1,2,8 | sed -e 's/ASM.*\+\;END=//g' > ager_mask.bed` 

Generate a "clean" vcf with: `perl nam_clean.pl ager.gvcf | grep "DP=27\|DP=25\|DP=26" | grep -v ":\.:" > ager_clean.vcf` (obviously changing file names so you're not writing over everything with `ager.vcf` etc. each time. You will want to change the depth filter to whatever you want. 

an play with the time windows (`-p` flag) if you want better resolution in some time frame, but the defaults are fine to start with.. You c
-bash: syntax error near unexpected token `('
-bash: cd: /group/jrigrp11/nam_gvcfs: No such file or directory
Command 'gatk' not found, did you mean:
  command 'gitk' from deb gitk (1:2.34.1-1ubuntu1.10)
  command 'gawk' from deb gawk (1:5.1.0-1ubuntu0.1)
Try: apt install <deb name>
ERROR:: command not found
-bash: syntax error near unexpected token `('
-bash: !/bin/bash: event not found
ERROR: Unable to locate a modulefile for 'GATK/4.2.3.0'
Autoloading htslib/1.13
Loading htslib/1.13

Loading samtools/1.13
  Loading requirement: htslib/1.13
-bash: cd: /group/jrigrp11/nam_gvcfs: No such file or directory
Command 'gatk' not found, did you mean:
  command 'gitk' from deb gitk (1:2.34.1-1ubuntu1.10)
  command 'gawk' from deb gawk (1:5.1.0-1ubuntu0.1)
Try: apt install <deb name>
-bash: command substitution: line 1: unexpected EOF while looking for matching `''
-bash: command substitution: line 2: syntax error: unexpected end of file
Generate: command not found
grep: ager.gvcf: No such file or directory
Command 'etc.' not found, did you mean:
  command 'etcd' from deb etcd-server (3.3.25+dfsg-7ubuntu0.22.04.1)
Try: apt install <deb name>
Generate: command not found
sed: can't read nam_clean.pl: No such file or directory
sed: can't read ager.gvcf: No such file or directory
**nam_clean.pl**: command not found
grep: warnings: No such file or directory
my: command not found
-bash: syntax error near unexpected token `NAM,'
my: command not found
my: command not found
-bash: syntax error near unexpected token `)'
++: command not found
-bash: syntax error near unexpected token `{'
-bash: syntax error near unexpected token `('
-bash: syntax error near unexpected token `('
Command 'next' not found, but can be installed with:
apt install mailutils-mh  # version 1:3.14-1, or
apt install mmh           # version 0.4-4
apt install nmh           # version 1.7.1-11
Ask your administrator to install one of them.
-bash: syntax error near unexpected token `('
-bash: syntax error near unexpected token `('
my: command not found
-bash: syntax error near unexpected token `('
^X^C
-bash: syntax error near unexpected token `}'
**nam_multihep.pl**: command not found
Command 'Combine' not found, did you mean:
  command 'combine' from deb moreutils (0.66-1)
Try: apt install <deb name>
-bash: command substitution: line 3: syntax error near unexpected token `('
-bash: command substitution: line 3: `Make msmc2 input file (format [here](https://github.com/stschiff/msmc-tools/blob/master/msmc-tutorial/guide.md)) '
-bash: command substitution: line 1: syntax error near unexpected token `then'
-bash: command substitution: line 1: ` then sort '
sort: cannot read: nam_multihep.pl: No such file or directory
gzip: ager_clean.vcf.gz: No such file or directory
gzip: ager_mask.bed.gz: No such file or directory
gzip: #.gz: No such file or directory
gzip: ager_all.txt.gz: No such file or directory
Command 'use' not found, did you mean:
  command 'muse' from deb muse (4.0.0-1build1)
  command 'ase' from deb ase (3.22.1-1ubuntu1)
  command 'nse' from deb ns2 (2.35+dfsg-3.1)
  command 'fuse' from deb fuse-emulator-gtk (1.6.0+dfsg1-2)
  command 'fuse' from deb fuse-emulator-sdl (1.6.0+dfsg1-2)
Try: apt install <deb name>
my: command not found
-bash: syntax error near unexpected token `FILE,'
my: command not found
my: command not found
-bash: syntax error near unexpected token `)'
Command 'chomp' not found, did you mean:
  command 'comp' from deb mailutils-mh (1:3.14-1)
  command 'comp' from deb mmh (0.4-4)
  command 'comp' from deb nmh (1.7.1-11)
Try: apt install <deb name>
Command 'next' not found, but can be installed with:
apt install mailutils-mh  # version 1:3.14-1, or
apt install mmh           # version 0.4-4
apt install nmh           # version 1.7.1-11
Ask your administrator to install one of them.
-bash: syntax error near unexpected token `('
my: command not found
my: command not found
my: command not found
-bash: syntax error near unexpected token `{'
-bash: syntax error near unexpected token `{'
my: command not found
my: command not found
-bash: syntax error near unexpected token `split'
Warning: unknown mime-type for "," -- using "application/octet-stream"
Warning: unknown mime-type for "\t," -- using "application/octet-stream"
Warning: unknown mime-type for "," -- using "application/octet-stream"
Warning: unknown mime-type for "\t," -- using "application/octet-stream"
Warning: unknown mime-type for "+1,\t," -- using "application/octet-stream"
Warning: unknown mime-type for "," -- using "application/octet-stream"
Warning: unknown mime-type for "," -- using "application/octet-stream"
Warning: unknown mime-type for "\n" -- using "application/octet-stream"
Error: no such file ","
Error: no such file "\t,"
Error: no such file ","
Error: no such file "\t,"
Error: no such file "+1,\t,"
Error: no such file ","
Error: no such file ","
Error: no such file "\n"
=0: command not found
-bash: syntax error near unexpected token `}'
else{: command not found
=++1-: command not found
-bash: syntax error near unexpected token `}'
=: command not found
-bash: syntax error near unexpected token `}'
> ^C
farm🌿:msmc2$ nano README.md
farm🌿:msmc2$ less README.md 
farm🌿:msmc2$ make genomic_import.sh
make: *** No rule to make target 'genomic_import.sh'.  Stop.
farm🌿:msmc2$ touch genomic_import.sh
farm🌿:msmc2$ less README.md 
farm🌿:msmc2$ nano genomic_import.sh 
farm🌿:msmc2$ less README.md 
farm🌿:msmc2$ touch select_variants.sh
farm🌿:msmc2$ less README.md 
farm🌿:msmc2$ nano #!/bin/bash
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
    -O nam_"$SLURM_ARRAY_TASK_ID
ERROR: Unable to locate a modulefile for 'GATK/4.2.3.0'
Loading samtools/1.13
  INFO: Module 'samtools/1.13' is already loaded
-bash: cd: /group/jrigrp11/nam_gvcfs: No such file or directory
> nano select_variants.sh 
> ^C
farm🌿:msmc2$ nano select_variants.sh 
farm🌿:msmc2$ less README.md 
farm🌿:msmc2$ touch nam_clean.pl
farm🌿:msmc2$ less README.md 

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

Generate file of sequenceable regions using `grep ";END=" ager.gvcf| cut -f 1,2,8 | sed -e 's/ASM.*\+\;END=//g' > ager_mask.bed` 

Generate a "clean" vcf with: `perl nam_clean.pl ager.gvcf | grep "DP=27\|DP=25\|DP=26" | grep -v ":\.:" > ager_clean.vcf` (obviously chang>

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
