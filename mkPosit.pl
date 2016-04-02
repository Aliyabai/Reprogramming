#!usr/bin/perl -w
use strict;
use Data::Dumper;

@ARGV == 1 or die "Usage: perl $0 <merge.txt>\n";
#my $symbol = $1 if ($ARGV[0] =~ /sample.(.*?).txt$/);


my %site;

open IN, "<$ARGV[0]" or die "$!";
<IN>;
while(<IN>){
        chomp;
        my @t = split /\t/;
        next if ($t[3] eq ".");
        my $base = (split /\-\>/,$t[3])[0];
        my $strand = ($t[1] eq $base)? "+":"-";
        my ($chr, $pos) = (split /_/, $t[0])[0,1];
        $site{$chr}{$pos} = $strand;
        $site{$chr}{$pos+1} = $strand;
        $site{$chr}{$pos-1} = $strand;
}
close IN;

open OUT, ">filter.pos.reg" or die "$!";
for my $x(1..22,"X","Y"){
        my $chr = "chr".$x;
        for my $pos(sort {$a <=> $b} keys %{$site{$chr}}){
                print OUT join("\t",$chr,$pos,$site{$chr}{$pos}),"\n";
        }
}
close OUT;


###merge.txt:
#Chr_Pos Ref     Alt     Mut_type        Genename        Strand  Alu     Region  SRR1107924|FHep SRR1107925|FHep SRR1107926|hiHep        SRR1107927|hiHep
 #       SRR1107928|hiHep        SRR1107931|HEF  SRR597900|HEF   SRR597904|HEF   SRR597908|HEF   SRR1107929|ES-Hep       SRR828795|ES-Hep        SRR828796|ES-Hep        SRR828797|ES-Hep        SRR828798|ES-Hep        SRR828799|ES-Hep        SRR828800|ES-Hep        SRR828801|ES-Hep        SRR828802|ES-Hep        SRR1107930|HepG2        SRR521470|HepG2 SRR521471|HepG2 SRR521472|HepG2 SRR521473|HepG2 SRR521474|HepG2 SRR521475|HepG2 SRR521476|HepG2
#chr10_101556205 A       G       A->G    ABCC2   +       Alu     intronic        0       0       0       0       0       0       0       0       0       0
#       0       0       0       0       0       0       0       0       0       0       0.272727272727273       0.230769230769231       0       0.166666666666667       0       0
