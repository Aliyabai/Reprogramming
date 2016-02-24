#!usr/bin/perl -w
use strict;
die "USAGE: perl $0 srrlist" if @ARGV < 1;
my $pileupStat = "bin/pileup.stat.pl";
my $dir = "/ifs5/ST_COMG/USER/baiyali/RNA_editing/Reprogramming/data_process/new.data.work/ADAR_motif";

open IN, "$ARGV[0]" or die $!;
open OUT, ">run2.sh" or die $!;
while(<IN>){
    chomp;
    my $srr = (split /\t/,)[0];
    print OUT "perl $pileupStat $dir/pipeup/$srr.pipeup > $dir/Stat/$srr.pipeup.stat\n";
}
close OUT;

###perl bin/pileup.stat.pl pipeup/SRR1107924.pipeup > Stat/SRR1107924.pipeup.stat
###建STAT文件
