#####Summary_1

#!usr/bin/perl -w
use strict;
die "usage: perl $0 SRRlist" if (@ARGV != 1);
my $maindir = "/ifs5/ST_COMG/USER/baiyali/RNA_editing/Reprogramming/data_process";
my $work = "/ifs5/ST_COMG/USER/baiyali/RNA_editing/Reprogramming/data_process/work/summary";


open IN,"$ARGV[0]" or die $!;
my @index = <IN>;
chomp @index;

for my $SRR (@index){
    open OUT, ">$work/$SRR.summary.sh" or die $!;
    my $fa = "$maindir/$SRR/clean/$SRR\_1.Clean.fq.gz";
    my $map_summary = "$maindir/$SRR/tophat/align_summary.txt";
    my $pre_dupfile = "$maindir/$SRR/bamprocess/$SRR\.reorder.bam";
    my $duppedfile = "$maindir/$SRR/bamprocess/$SRR\.reorder.markdup.bam";
    my $final_bam = "$maindir/$SRR/bamprocess/$SRR\.reorder.markdup.baseRecal.bam";
    print OUT "echo fa_reads:;less $fa |wc -l\n";
    print OUT "echo mapped_reads:;grep Mapped $map_summary |awk '{print \$3}'\n ";
    print OUT "echo pre_dupped:;samtools view -h $pre_dupfile |grep SRR|grep -v ID:|wc -l \n";
    print OUT "echo dupped:;samtools view $duppedfile |grep SRR|grep -v ID:|wc -l\n";
    print OUT "echo uniqmap_reads:;samtools view $final_bam |grep NH:i:1 |grep SRR |grep -v ID:|wc -l";
    close OUT;
}

##### Summary_2 #####
#!usr/bin/perl -w 
use strict;
die "USAGE:perl $0 SRRlist" if (@ARGV != 1);
my $miandir = "/ifs5/ST_COMG/USER/baiyali/RNA_editing/Reprogramming/data_process";
my $work = "$miandir/work/summary";

open IN,"$ARGV[0]" or die $!;
my @index = <IN>;
chomp @index;

open OUT, ">$work/summary.txt" or die $!;
print OUT "SRR\tPE/SE\tRead_length\tClean_reads\tClean_bases\tMapped_reads\tDuplicated_reads\tDuplication_rates\tUniq_mapped_reads\tUniq_mapped_rates\n";
for my $SRR (@index){
    my @summary = `less $work/$SRR.summary.sh.o* `;
    chomp @summary;
    my $fa_reads = $summary[1];
    my $mapped_reads = $summary[3];
    my $pre_dupped = $summary[5];
    my $dupped = $summary[7];
    my $uniq = $summary[9];
    
    my $clean_reads = $fa_reads/4;
    my $clean_bases = $clean_reads * 100;
    my $Duplicated_reads = $pre_dupped - $dupped;
    my $Duplication_rates = $Duplicated_reads/$pre_dupped;
    my $Uniq_mapped_rates = $uniq / $dupped;
    print OUT "$SRR\tSE\t100bp\t$clean_reads\t$clean_bases\t$mapped_reads\t$Duplicated_reads\t$Uniq_mapped_rates\t$uniq\t$Uniq_mapped_rates\n";

}
close OUT;

