#!usr/bin/perl -w
use strict;
die "perl $0 srrlist" if (@ARGV != 1);
my $maindir = "/ifs5/ST_COMG/USER/baiyali/RNA_editing/Reprogramming/data_process";
my $fadir = "/ifs5/ST_COMG/USER/baiyali/RNA_editing/Reprogramming/raw_data/fa";
open IN, "$ARGV[0]" or die $!;
my @srrlist = <IN>;
chomp @srrlist;

for my $srr (@srrlist){
    unless (-d "$maindir/$srr"){`mkdir $maindir/$srr `;}
    unless (-d "$maindir/$srr/clean"){`mkdir $maindir/$srr/clean `;}
    unless (-d "$maindir/$srr/tophat"){`mkdir $maindir/$srr/tophat `;}
    open OUT, ">/ifs5/ST_COMG/USER/baiyali/RNA_editing/Reprogramming/data_process/work/clean_tophat/$srr.tophat_only.sh"or die $!;
    print OUT "/ifs5/ST_ANNO/USER/chenrongchang/Tools/SOAPnuke1.5.0 filter -1 $fadir/$srr\_1.fastq.gz -o $maindir/$srr/clean -C $maindir/$srr/clean/$srr\_1.Clean.fq.gz;echo 1 complete\n";
    print OUT "export PATH=/ifs1/ST_SINGLECELL/USER/jiangrunze/tool/tophat-2.0.12.Linux_x86_64/bin/:\$PATH\n";
    print OUT "/ifs1/ST_SYSBIO/USER/liudb2/RNA_editing/bladder/xiongheng/bin/tophat-2.0.12.Linux_x86_64/tophat2 --solexa-quals --read-mismatches 2 --read-gap-length 3 --read-edit-dist 3 --library-type fr-unstranded -p 6 -r 30 --b2-fast --rg-center bgi --rg-platform illumina --no-novel-juncs --no-novel-indels --transcriptome-index=/ifs1/ST_SYSBIO/USER/liudb2/RNA_editing/bladder/xiongheng/database/transcriptomeindex/tranindex  --rg-id $srr --rg-sample $srr --rg-library D276HACXX_L5 -o $maindir/$srr/tophat /ifs1/ST_SYSBIO/USER/liudb2/RNA_editing/bladder/xiongheng/database/hg19index/hg19 $maindir/$srr/clean/$srr\_1.Clean.fq.gz;echo 2 complete";
    close OUT;
}
