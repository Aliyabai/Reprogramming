#!usr/bin/perl -w
use strict;
die "perl $0 srrlist" if (@ARGV != 1);
my $maindir = "/ifs5/ST_COMG/USER/baiyali/RNA_editing/Reprogramming/raw_data/work";
my $sradir = "/ifs5/ST_COMG/USER/baiyali/RNA_editing/Reprogramming/raw_data/SRR";
my $fastq_dump = "/ifs1/ST_CANCER/USER/zhouyanqing/bin/sratoolkit.2.4.1-ubuntu64/bin/fastq-dump";
open IN, "$ARGV[0]" or die $!;
my @srrlist = <IN>;
chomp @srrlist;

for my $srr (@srrlist){
    my $sra = "$sradir/$srr/$srr.sra";
    my $fa = "/ifs5/ST_COMG/USER/baiyali/RNA_editing/Reprogramming/raw_data/fa";
    open OUT, ">$maindir/$srr.clean_tophat.sh"or die $!;
    print OUT "$fastq_dump $sra --split-files -O $fa;echo 1 complete\n";
    print OUT "perl /ifs5/ST_COMG/USER/liwenhui/RNA_editing/ADARs/pan-cancer/pan/3.SingleCell/bin/id.modf1.pl $fa/$srr\_1.fastq $fa/$srr\_1.fq\nrm $fa/$srr\_1.fastq;echo 2 complete\n";
    print OUT "perl /ifs5/ST_COMG/USER/liwenhui/RNA_editing/ADARs/pan-cancer/pan/3.SingleCell/bin/id.modf2.pl $fa/$srr\_2.fastq $fa/$srr\_2.fq\nrm $fa/$srr\_2.fastq;echo 3 complete";
    close OUT;
}
