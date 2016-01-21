use strict;
use warnings;

@ARGV == 1 or die "perl $0 <pileup> > output\n";

if ($ARGV[0] =~ /\.gz$/){
        open IN,"gzip -dc $ARGV[0] |" or die $!;
}else{
        open IN,"$ARGV[0]" or die $!;
}

print"#Chr\tPos\tDP\tRef\tRef_Sup\tAlt\tAlt_Sup\n";
while(<IN>){
        chomp;
        my @a = split; #$chr pos ref effective_reads bases quals mapQs tposs misms qlens
        my $ref = uc($a[2]);
        my ($dp,$ref_sup,$base,$alt_sup) = ($a[3],0,"N",0);
        if($a[4]){
                $base =uc($a[4]);
                $ref_sup += ($base =~ s/$ref//g);
                $alt_sup = length($base);
                if($base){
                        my %base;
                        foreach(split //,$base){
                                $base{$_}++;
                        }
                        $base = (sort {$base{$b} <=> $base{$a}}keys %base)[0];
                        $alt_sup = $base{$base};
                        $dp = $ref_sup + $alt_sup;
                }else{
                        $base = "N";
                }
        }else{
                print STDERR "Warning:$_\n" if($a[3] != 0);
        }       
        print"$a[0]\t$a[1]\t$dp\t$ref\t$ref_sup\t$base\t$alt_sup\n";
}
close IN;
