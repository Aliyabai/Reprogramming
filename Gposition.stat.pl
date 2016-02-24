#!usr/bin/perl -w
use strict;
use Data::Dumper;

@ARGV == 4 or die "Usage: perl $0 <PipeupSite.stat.list> <pos.reg> <merge.txt> <Base>\n";
my $base = pop @ARGV;
my %list = qw{G C C G A T T A};

my $antibase = $list{$base};
my $type = $1 if ($ARGV[2] =~ /merge.(.*?).sample.txt$/); print "$type\n";
%strand = &Strand($ARGV[1]);
%site = &Site($ARGV[2]);

my %data2;
for my $sid(sort keys %pfile){ #print "$sid\n";
        my %reg = &pipeupStat($sid);
        for my $x(1..22,"X","Y"){
                my $chr = "chr".$x;
                next unless (exists $site{$sid}{$chr});
                for my $pos(sort {$a <=> $b} keys %{$site{$sid}{$chr}}){
                        my $strand = $site{$sid}{$chr}{$pos};
                        if ($strand eq "+"){
#                               $data{$sid}{upstream}++ if ($reg{$chr}{$pos-1} eq $base);
#                               $data{$sid}{downstream}++ if ($reg{$chr}{$pos+1} eq $base);
                                $data2{$sid}{unstream}{$reg{$chr}{$pos-1}}++ if (exists $reg{$chr}{$pos-1});
                                $data2{$sid}{downstream}{$reg{$chr}{$pos+1}}++ if (exists $reg{$chr}{$pos+1});
                        }else{
#                               $data{$sid}{upstream}++ if ($reg{$chr}{$pos+1} eq $antibase);
#                               $data{$sid}{downstream}++ if ($reg{$chr}{$pos-1} eq $antibase);
                                $data2{$sid}{unstream}{$reg{$chr}{$pos+1}}++ if (exists $reg{$chr}{$pos+1});
                                $data2{$sid}{downstream}{$reg{$chr}{$pos-1}}++ if (exists $reg{$chr}{$pos-1});
                        }
                        $data{$sid}{total}++;
                        $data2{$sid}{total}++;
                }#print "$data2{$sid}{total}\n";
        }
}

my @bases = qw{A T C G};
open OUT, ">out/$type.percentage.stat" or die "$!";
for my $sid(@order){
        if (!exists $data2{$sid}){
                print OUT $sid,"\n";
                next;
        }
        my $total = $data2{$sid}{total};
        my @line = ($sid);
        for my $region("unstream","downstream"){
                for my $b(@bases){
                        my $tmp = (exists $data2{$sid}{$region}{$b})? $data2{$sid}{$region}{$b} : 0;
                        push @line, ($tmp,$tmp/$total);
                }
        }
        print OUT join("\t",@line),"\n";
}
close OUT;

##############################################################################
sub pipeupStat{
        my $ssid = shift; #print "$ssid\n";
        my $input = $pfile{$ssid}; #print "$input\n";
        my %hash;
        open FI, "<$input" or die "$!";
        while(my $line = <FI>){
                next if ($line =~ /^#Chr/);
                chomp $line;
                my @t = split /\t/,$line;
                next unless (exists $site{$ssid}{$t[0]}{$t[1]} or exists $site{$ssid}{$t[0]}{$t[1]+1} or exists $site{$ssid}{$t[0]}{$t[1]-1});
                next if ($t[2] == 0);
                my $strand = $strand{$t[0]}{$t[1]};
                if ($t[5] eq "N"){
                        if ($strand eq "+"){
                                $hash{$t[0]}{$t[1]} = $t[3];
                        }else{
                                $hash{$t[0]}{$t[1]} = $list{$t[3]};
                        }
                }else{
                        if ($t[4] == $t[6]){
                                if ($strand eq "+"){
                                        if ($t[3] eq $base or $t[5] eq $base){
                                                $hash{$t[0]}{$t[1]} = $base;
                                        }else{
                                                $hash{$t[0]}{$t[1]} = $t[3];
                                        }
                                }else{
                                        if ($t[3] eq $antibase or $t[5] eq $antibase){
                                                $hash{$t[0]}{$t[1]} = $base;
                                        }else{
                                                $hash{$t[0]}{$t[1]} = $list{$t[3]};
                                        }
                                }
                        }else{
                                my $current_base = ($t[4]>$t[6])? $t[3]:$t[5];
                                if ($strand eq "+"){
                                        $hash{$t[0]}{$t[1]} = $current_base;
                                }else{
                                        $hash{$t[0]}{$t[1]} = $list{$current_base};
                                }
                        }
                }
        }
        close FI;
        return %hash;
}

sub Site{
        my $input = shift;
        my %hash;
        open FI, "<$input" or die "$!";
        my $head = <FI>; chomp $head; my @head = split /\t/,$head; #print "@head\n";
        while(my $line = <FI>){
                chomp $line;
                my @t = split /\t/,$line;
                next if ($t[3] eq ".");
                my ($chr, $pos) = (split /_/,$t[0])[0,1];
                my $strand = $strand{$chr}{$pos};
                for my $i(8..$#t){                                    # print "$i\n";
                        next if ($t[$i] == 0);#print "$t[$i]\n";
                        my $sid = (split /\|/,$head[$i])[0]; #print "$sid\n";
                        $hash{$sid}{$chr}{$pos} = $strand; #print "$strand\t$hash{$sid}{$chr}{$pos}\n";
                }
        }
        close FI;
        return %hash;
}

sub Strand{
        my $input = shift;
        my %hash;
        open FI, "<$input" or die "$!";
        while(my $line = <FI>){
                chomp $line;
                my ($chr,$pos,$strand) = (split /\t/,$line)[0..2];
                $hash{$chr}{$pos} = $strand; #print "$chr\t$pos\t$strand\t$hash{$chr}{$pos}\n";
        }
        close FI;
        return %hash;
}

sub read_file{
        my $input = shift;
        my %hash;
        open FI, "<$input" or die "$!";
        while(my $line = <FI>){
                chomp $line;
                my ($sid,$file) = (split /\t/,$line)[0,1];
                $hash{$sid} = $file;  #print "$sid\t$hash{$sid}\n";
        }
        close FI;
        return %hash;
}
