#!usr/bin/perl -w
use strict;
die "perl $0 <in.bam> <out1.fq> <out2.fq> <outs.fq> <region.txt> <vreadsnum>" unless @ARGV==6;
my $input=shift @ARGV;
my $out1=shift @ARGV;
my $out2 = shift @ARGV;
my $outs = shift @ARGV;
my $inregion = shift @ARGV;
my $numout = shift @ARGV;
open REGION, $inregion or die $!;
open OUT, ">$numout"or die $!;
open OUT1, ">$out1" or die $!;
open OUT2, ">$out2" or die $!;
open OUTS, ">$outs" or die $!;
my %hash;
while (<REGION>){
    chomp;
    my @array = split /\t/,$_;
    my $chr = $array[0];
    if ($chr eq "chrM_rCRS") {next;}
    my $refPos = $array[1];
    my $refbase = substr($array[3],0,1);
    my $samPos = "$chr:$refPos-$refPos";
    my $altbase = $array[11];
    &reads ($refPos,$samPos,$refbase,$altbase)
}
sub reads {
    my ($refPos,$samPos,$refbase,$altbase) = @_;
    open IN, "samtools view  $input $samPos |" or die $!;
    my $i = 0;
    while (<IN>){
        chomp;
        my @lines = split /\s+/,$_;
        my $dup = 1;
        my $uniq = 2;
        if ($lines[1] & 0x0400){$dup = 2;}
        if ($_ =~ /X0:i:1\t/){$uniq = 1;}
        if ($_ =~ /NH:i:1\t/){$uniq = 1;}
        if ($_ !~ /^@/ and  $dup ==1 and $uniq ==1  and $lines[4]>=20){
            my $name = $lines[0];
            if ($lines[0] =~ /\/1$/ || $lines[0] =~ /\/2$/){
                my $len = length($lines[0])-2;
                $name = substr($lines[0],0,$len);
            }
            my $flag = $lines[1];
            my $chr = $lines[2];
            my $alignPos = $lines[3];
            my @str = split /\d+/,$lines[5];
            shift @str;
            my @num = split /[A-Z]/,$lines[5];
            my $length = @str;
            my @cigar;
            for (my $i=0;$i<$length;$i++){
                my @subarray = ($str[$i],$num[$i]);
                push @cigar, [@subarray];
            }
            my $pchr = $lines[6];
            my $palignPos = $lines[7];
            my $seq = $lines[9];
            my $qual = $lines[10];
            my $end;
            my $readPos;
            my $signal = " ";
            for (my $c=0; $c<@cigar; $c++) {
                my $char=$cigar[$c][0];
                my $cnum=$cigar[$c][1];
                if($char eq "S"){
                    $readPos+=$cnum;
                }elsif($char eq "H"){
                    next;
                }elsif($char eq "M"){
                    $end=$alignPos+$cnum-1;
                    if($end<$refPos){
                        $readPos+=$cnum;
                        $alignPos+=$cnum;
                    }elsif($end>=$refPos){
                        $readPos+=$refPos-$alignPos+1;
                        last;
                    }
                }elsif($char eq "I"){
                    $readPos+=$cnum;
                }elsif($char eq "D"){
                    $end=$alignPos+$cnum-1;
                    if($end<$refPos){
                        $alignPos+=$cnum;
                    }else{
                         last;
                    }
                }elsif($char eq "N"){
                     $end=$alignPos+$cnum-1;
                     if($end<$refPos){
                         $alignPos+=$cnum;
                     }else{
                           $signal = "down";
                           last;
                      }
                }else{
                       print STDERR "Cigar process error at $alignPos\n";
                     }
            }
            if ($signal ne "down"){
                my $qbase  = substr($seq,$readPos-1,1);
                my $basequa = substr($qual,$readPos-1,1);
                my $basequat = ord($basequa)-33;
                if ($qbase ne $refbase and $altbase eq $qbase and $basequat >=20){
                    my $bin=unpack("B32",pack("N",$flag));
                    my @bin=split //,$bin;
                    @bin=reverse @bin;
                    if($bin[4]==1){
                        my @seq=split //,$seq;
                        my @qual=split //,$qual;
                        @seq=reverse @seq;
                        @qual=reverse @qual;
                        for (my $i=0; $i<@seq ; $i++){
                            $seq[$i]=~tr/ACGTacgt/TGCAtgca/;
                        }
                        $seq=join '',@seq;
                        $qual=join '',@qual;
                        $readPos = length($seq)+1-$readPos;
                     }
                     if ($flag & 0x0008 or $pchr eq "\*"){
                         print OUTS "\@".$name."_".$samPos."-".$readPos."\n".$seq."\n"."+"."\n".$qual."\n";
                         $i += 1;
                     }
                     else {
                         my $pe;
                         if ($bin[6]==1){$pe = 1;}
                         if ($bin[7]==1){$pe = 2;}
                         my $pep;
                         if ($bin[6]==1){$pep = 2;}
                         if ($bin[7]==1){$pep = 1;}
                         $i += 1;
                         $hash{$name}{$pep}{$samPos} = "\@".$name."_".$samPos."-\t$readPos\t-".$pe."\t".$seq."\n+\n".$qual."\n";
                     }
                }
            }
        }
    }
    my $outpos = (split /-/,$samPos)[1];
    my $outchr = (split /:/,$samPos)[0];
    print OUT $outchr."\t".$outpos."\t".$i."\n";
}
close IN;
my @keys = keys %hash;
if (@keys > 0){
    open IN1, "samtools view  $input|" or die $!;
    while (<IN1>){
        chomp;
        my @lines = split /\s+/,$_;
        my $dup = 1;
        my $uniq = 2;
        if ($_ =~ /X0:i:1\t/){$uniq = 1;}
        if ($_ =~ /NH:i:1\t/){$uniq = 1;}
        my $name = $lines[0];
        if ($lines[0] =~ /\/1$/ || $lines[0] =~ /\/2$/){
            my $len = length($lines[0])-2;
            $name = substr($lines[0],0,$len);
        }
        if ($lines[1] & 0x0400){$dup = 2;}
        if ($_ !~ /^@/ and $dup ==1 and $uniq == 1  and $lines[4]>=20){
            my $bin=unpack("B32",pack("N",$lines[1]));
            my @bin=split //,$bin;
            @bin=reverse @bin;
            my $pep;
            if ($bin[6]==1){$pep = 1;}
            if ($bin[7]==1){$pep = 2;}
            if (exists $hash{$name}{$pep}){
                if($bin[4]==1){
                    my @seq=split //,$lines[9];
                    my @qual=split //,$lines[10];
                    @seq=reverse @seq;
                    @qual=reverse @qual;
                    for (my $i=0; $i<@seq ; $i++){
                        $seq[$i]=~tr/ACGTacgt/TGCAtgca/;
                    }
                    $lines[9]=join '',@seq;
                    $lines[10]=join '',@qual;
                }
                for my $key (keys %{$hash{$name}{$pep}}){
                    my @array = split /\t/,$hash{$name}{$pep}{$key};
                    my $pe = $array[2];
                    my $readPos = $array[1];
                    if ($pep == 1){
                        print OUT1 "\@".$name."_".$key."-".$readPos.$pe."-1"."\n".$lines[9]."\n"."+"."\n".$lines[10]."\n";
                        print OUT2 $array[0].$array[1].$array[2]."-2"."\n".$array[3];
                    }
                    if ($pep == 2){
                        print OUT2 "\@".$name."_".$key."-".$readPos.$pe."-2"."\n".$lines[9]."\n"."+"."\n".$lines[10]."\n";
                        print OUT1 $array[0].$array[1].$array[2]."-1"."\n".$array[3];
                    }
                }
            }
        }
    }
}
close OUT1;
close OUT2;
close OUTS;
close IN1;
