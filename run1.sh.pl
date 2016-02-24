#!usr/bin/perl -w
use strict;
die "USAGE: perl $0 srrlist" if @ARGV < 1;
my $pileupSite = "/bin/pileupSite.pl";
my $ref = "database/reference/hg19/hg19.fa";
my $reg = "filter.pos.reg";

open IN, "$ARGV[0]" or die $!;
open OUT, ">run1.sh" or die $!;
while(<IN>){
    chomp;
    my ($srr, $maindir) = (split /\t/,)[0,2];
    my $bam = "$maindir/$srr.reorder.markdup.baseRecal.bam";
    print OUT "/usr/bin/perl $pileupSite --bam $bam --reg $reg --out /ifs5/ST_COMG/USER/baiyali/RNA_editing/Reprogramming/data_process/new.data.work/ADAR_motif/pipeup/$srr.pipeup --ref $ref --minMQ 20 --minBQ 20 --uniq\n";
}
close OUT;
