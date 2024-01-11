use strict;
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
