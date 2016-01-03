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
        next if ($t[6] eq ".");
        my $base = (split /\-\>/,$t[6])[0];
        my $strand = ($t[4] eq $base)? "+":"-";
        $site{$t[0]}{$t[1]} = $strand;
        $site{$t[0]}{$t[1]+1} = $strand;
        $site{$t[0]}{$t[1]-1} = $strand;
}
close IN;

open OUT, ">Trim15bp.filter2.pos.reg" or die "$!";
for my $x(1..22,"X","Y"){
        my $chr = "chr".$x;
        for my $pos(sort {$a <=> $b} keys %{$site{$chr}}){
                print OUT join("\t",$chr,$pos,$site{$chr}{$pos}),"\n";
        }
}
close OUT;
