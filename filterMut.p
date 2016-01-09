#!usr/bin/perl -w
use strict;
die "perl $0 <in.mut.txt.gz> \n"if(@ARGV!=1);
open IN,"gzip -dc $ARGV[0]|" or die $!;
while(<IN>){
    chomp;
    if($_=~/^Chr/){next;}
    my @c=split /\s+/,$_;
    my $freq=$c[12]/$c[2];
    my $sb=$c[15]/($c[15]+$c[16]);
    if($c[12]>=$c[19] && $freq>=0.1 && $sb>=0.1 && $sb<=0.9){
        print "$_\n";
    }
}
close IN;


#输入文件格式：
#Chr     Pos     Depth   Ref     Sup     Mid     End     Plus    Minus   Mism    NoMism  Alt     Sup     Mid     End     Plus    Minus   Mism    NoMism     rposPvalue      sbPvalue        mmPvalue
#chr1    14677   11      G       4       3       1       2       2       1       3       A       7       6       1       3       4       2       5       2      1       1       1

