#!usr/bin/perl -w
use strict;
my ($in,$vreadsnum,$out) = @ARGV;
open IN, $in or die $!;
open OUT, ">$out" or die $!;
while (<IN>){
    chomp;
    my @array = split /\t/,$_;
    if ($array[9] eq "Alu" && $array[6] > $vreadsnum){
        print OUT $_."\n";
    }
    if ($array[9] eq "Non-Alu" && $array[6] > $vreadsnum && $_ !~ /splicing/){
        print OUT $_."\n";
    }
}
close OUT;
close IN;
