#!/usr/bin/perl -w
use strict;
my ($in,$out) = @ARGV;
open IN, $in or die $!;
open OUT, ">$out" or die $!;
while (<IN>){
    chomp;
    if ($_ =~ /^\@/ and $_ =~ /length=100$/){
        my @a = split /\s+/,$_;
        pop(@a);
        my $line = join ("_",@a);
        print OUT $line."/1\n";
    }
    elsif ($_ =~ /^\+/ and $_ =~ /length=100$/){
        my @b = split //,$_;
        print OUT $b[0]."\n";
    }
    else {print OUT $_."\n";}
}
close OUT;
close IN;
