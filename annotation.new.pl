#!usr/bin/perl -w
use strict;
use FindBin qw($Bin);
die "perl $0 <in.mutDet.result> <alu.bed> <outprefix>\n"if(@ARGV!=3);
my ($in,$aluBed,$prefix)=@ARGV;

die "$Bin/annovar/annotate_variation.pl does not exits!\n" if(!-f "$Bin/annovar/annotate_variation.pl");
die "$Bin/annovar/humandb_hg19 does not exists!\n" if(!-d "$Bin/annovar/humandb_hg19");
my %hash=();
my $n=0;
my $NorT="Normal";
if($in=~/Tumor/i){$NorT="Tumor";}
print "Reading input file\n";
open IN,"$in" or die $!;
open A,">$prefix.anno" or die $!;
while(<IN>){
    chomp;
    if($_=~/#/){next;}
    my @c=split /\s+/,$_;
    @{$hash{$c[0]}{$c[1]}}=($c[3],$c[11],".",".",".","Non-Alu",$c[2],$c[4],$c[12]); # ref,alt,strand,region,gene,alu,depth,ref-sup,alt-sup
    print A "$c[0]\t$c[1]\t$c[1]\t$c[3]\t$c[11]\n";
    $n++;
}
close IN;
close A;
print "$n sites are loaded\n";

`perl $Bin/annovar/annotate_variation.pl --splicing_threshold 5 -buildver hg19 $prefix.anno $Bin/annovar/humandb_hg19`;

$n=0;
print "Reading annotation file\n";
open ANNO,"$prefix.anno.variant_function" or die $!;
while(<ANNO>){
    chomp;
    my @c=split /\s+/,$_;
    my $ch=$c[2];
    my $ps=$c[3];
    if(exists $hash{$ch}{$ps}){
        $hash{$ch}{$ps}[3]=$c[0];
        if($c[0]!~/intergenic/){
            my $g="";
            if($c[0] eq "splicing" || $c[0] eq "ncRNA_splicing"){
                if($c[1]=~/\(/){
                    $g=(split /\(/,$c[1])[0];
                }else{
                    $g=$c[1];
                }
            }else{
                $g=$c[1];
            }
            if($g=~/,/){
                $g=(split /,/,$g)[0];
            }
            if($g=~/;/){
                $g=(split /;/,$g)[0];
            }
            $hash{$ch}{$ps}[4]=$g;
            $n++;
        }else{
            $hash{$ch}{$ps}[4]=$c[1];
        }
    }
}
close ANNO;
print "$n sites are in non intergenic region\n";

$n=0;
print "Reading refGene file\n";
my %gene=();
open RG,"$Bin/annovar/humandb_hg19/hg19_refGene.txt" or die $!;
while(<RG>){
    chomp;
    my @c=split /\s+/,$_;
    $gene{$c[-4]}=$c[3];
    $n++;
}
close RG;
print "$n transcripts are loaded\n";

$n=0;
print "Alu region annotation\n";
`perl $Bin/pos_filterByRegion.pl $aluBed $in >$prefix.alu.txt`;
open ALU,"$prefix.alu.txt" or die $!;
while(<ALU>){
    chomp;
    if($_=~/#/){next;}
    my @c=split /\s+/,$_;
    if(exists $hash{$c[0]}{$c[1]}){
        $hash{$c[0]}{$c[1]}[5]="Alu";
        $n++;
    }
}
close ALU;
print "$n sites are in Alu regions\n";

$n=0;
print "Out put results\n";
open OUT,">$prefix.anno.txt" or die $!;
foreach my $chr(sort keys %hash){
    foreach my $pos(sort {$a<=>$b} keys %{$hash{$chr}}){
        my $ref=$hash{$chr}{$pos}[0];
        my $alt=$hash{$chr}{$pos}[1];
        if($hash{$chr}{$pos}[3] =~ /intergenic/){
            print OUT "$chr\t$pos\t$ref\t$alt\t$hash{$chr}{$pos}[6]\t$hash{$chr}{$pos}[7]\t$hash{$chr}{$pos}[8]\t.\t.\t$hash{$chr}{$pos}[5]\t$ha
sh{$chr}{$pos}[4]\t$hash{$chr}{$pos}[3]\n";
        }else{
            my $gene=$hash{$chr}{$pos}[4];
            my $strand=$hash{$chr}{$pos}[2];
            if($hash{$chr}{$pos}[4]=~/(\w+),/){
                $gene=$1;
            }
            if(exists $gene{$gene}){
                $strand=$gene{$gene};
            }else{
                print "Gene $gene does not exists in refGene\n";
            }
            $ref=~tr/ATCG/TAGC/ if( $strand eq "-" );
            $alt =~ tr/ATCG/TAGC/ if( $strand eq "-" );
            my $type = "$ref->$alt";
            print OUT "$chr\t$pos\t$ref\t$alt\t$hash{$chr}{$pos}[6]\t$hash{$chr}{$pos}[7]\t$hash{$chr}{$pos}[8]\t$type\t$strand\t$hash{$chr}{$po
s}[5]\t$gene\t$hash{$chr}{$pos}[3]\n";
            $n++;
        }
    }
}
close OUT;
print "$n sites are not in intergenic regions\n";

`cut -f8 $prefix.anno.txt |sort|uniq -c >$prefix.stat.txt`;
#`head -n 12 $prefix.stat.txt|awk '{a+=\$1}END{print a}' >>$prefix.stat.txt`;
`cut -f10 $prefix.anno.txt |sort|uniq -c >>$prefix.stat.txt`;
`cut -f12 $prefix.anno.txt |sort|uniq -c >>$prefix.stat.txt`;
`perl $Bin/count.pl $prefix.stat.txt >$prefix.stat.xls`;
#`rm -f $prefix.anno $prefix.anno.log $prefix.anno.exonic_variant_function $prefix.anno.variant_function $prefix.alu.txt`;
