#!usr/bin/perl -w
use strict;
die "USAGE: perl $0 srrlist" if @ARGV <1;
open IN, "$ARGV[0]" or die $!;
my @SRR = <IN>;
chomp @SRR;
my $main = "/ifs5/ST_COMG/USER/baiyali/RNA_editing/Reprogramming/data_process";

foreach my $srr(@SRR){
   my $bam = "$main/$srr/bamprocess/$srr.reorder.markdup.baseRecal.bam";   ###bbamamprocess 处理后的文件
   open OUT, ">$main/work/RPKM/R/$srr.R" or die "$!";
   print OUT "library(methods)\nlibrary(Rsubread)\nlibrary(limma)\nlibrary(edgeR)\nfc <- featureCounts(files=\"$bam\",isPairedEnd=T,annot.inbuilt=\"hg19\")\nwrite.table(fc\$counts,file=\"$main/work/RPKM/R/$srr.readscount\",sep=\"\\t\")\nx <- DGEList(counts=fc\$counts,genes=fc\$annotation[,c(\"GeneID\",\"Length\")])\nx_rpkm <- rpkm(x,x\$genes\$Length)\nwrite.table(x_rpkm,file=\"$main/work/RPKM/R/$srr..rpkm.txt\",sep=\"\\t\")\n";   
   close OUT;
   open OUT, ">$main/work/RPKM/sh/$srr.sh" or die $!;
   print OUT "/ifs1/ST_SINGLECELL/USER/yuqichao/tools/R-3.1.3/bin/Rscript\n$main/work/RPKM/R/$srr.R ; echo 1 is OK! ";
   close OUT;
}
