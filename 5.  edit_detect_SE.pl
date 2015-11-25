#!usr/bin/perl -w      
use strict;
die "perl $0 rnabamlist outdir" if (@ARGV != 2);
my ($rnabamlist,$outdir) = @ARGV;
my $bin = "/ifs1/ST_SYSBIO/USER/liudb2/RNA_editing/bladder/xiongheng/bin";
my $chr = "/ifs1/ST_SYSBIO/USER/liudb2/database/reference/hg19/hg19.fa";
my $refgene = "/ifs1/ST_SYSBIO/USER/liudb2/database/annovar/humandb_hg19/hg19_refGene.txt";
my $database = "/ifs1/ST_SYSBIO/USER/liudb2/RNA_editing/bladder/xiongheng/database/datanew/hg_mrna.fa";
my $bwa = "/ifs1/ST_SYSBIO/USER/liudb2/bin/TumorCare/bin/bwa";
my $dbsnp = "/ifs1/ST_SYSBIO/USER/liudb2/database/GATK/2.8/hg19/dbsnp_138.hg19.vcf.gz";
my $simpleRepeat = "/ifs1/ST_SYSBIO/USER/liudb2/bin/TumorCare/database/repeat/hg19_simpleRepeat.reg.bed";
my $snpdir = "/ifs1/ST_SYSBIO/USER/liudb2/RNA_editing/bladder/xiongheng/dbsnp";
open IN, $rnabamlist or die $!;
my @list;
my %hash;
while (<IN>){
    chomp;
    my @array = split /\s+/,$_;
    $hash{$array[0]} = $array[-1];
}
close IN;
foreach my $SRR (sort keys %hash){
   
        my $bam = $hash{$SRR};
        unless (-d "$outdir/$SRR/editing_detect") {`mkdir -p $outdir/$SRR/editing_detect`;}
        my $output = "$outdir/$SRR/editing_detect";
        open OUT, ">$outdir/work/edit_detect/$SRR\.removesnp2.sh" or die $!;
        ##MismatchStat##mismatchrate -x using the first int reads to calculate      -u use unique mapping reads to calculate mismatch rate
        print OUT "$bin/MismatchStat -i $bam -x 1000000 -u -o $output/stat.txt;echo 1 complete\n";
        ##MutDet## -r [str] reference fasta file   -q [int] minimal mapping quality   -v [float] file including the backgroud mismatch rate    -u disc
ard the non-unique mapping reads 筛选最小匹配率大于20包含mismatch率的unique reads
        print OUT "$bin/MutDet -i $bam -r $chr -q 20 -v $output/stat.txt -u -o $output/binomFilter.txt.gz;echo 2 complete\n";
        #$freq=$c[12]sup/$c[2]depth;$sb=$c[15]plus/($c[15]plus+$c[16]minus);if $c[12]>=$c[19]index && $freq>=0.1 && $sb>=0.1 && $sb<=0.9# print
        print OUT "perl $bin/filterMut.pl $output/binomFilter.txt.gz > $output/binomFilter.filt.fre.strand.txt;echo 3 complete\n";
        #
        print OUT "perl $bin/end.strand.mism.poly.rep.filt.pl --in $output/binomFilter.filt.fre.strand.txt --out $output/binomFilter.filt.fre.strand.
end.simplerepeat.homopolymer.txt --simpleRepeat $simpleRepeat --reference $chr --endfre 0.9 --endp 0.05 --strandfre 0.9 --strandp 0.005 --homopolymer 
5;echo 4 complete\n";
        #
       print OUT "perl $bin/pos_filter.pl --input $output/binomFilter.filt.fre.strand.end.simplerepeat.homopolymer.txt --loci $dbsnp --filt 2 --out 
$output/binomFilter.filt.fre.strand.end.simplerepeat.homopolymer.dbsnp.txt;echo 5 complete\n";
        #
        print OUT "perl $bin/pos_filter.pl --input $snpdir/1000g.YH.Mongolian.dai.Dutch.snp.txt.gz --loci $output/binomFilter.filt.fre.strand.end.sim
plerepeat.homopolymer.dbsnp.txt --out $output/1000g.YH.Mongolian.dai.Dutch.exchange.txt;echo 6 complete\n";
        #
        print OUT "perl $bin/pos_filter.pl --input $output/binomFilter.filt.fre.strand.end.simplerepeat.homopolymer.dbsnp.txt -loci $output/1000g.YH.
Mongolian.dai.Dutch.exchange.txt --filt 2 --out $output/binomFilter.filt.fre.strand.end.simplerepeat.homopolymer.dbsnp.snp.txt;echo complete 7\n";


        unless (-d "$output/bwa") {`mkdir -p $output/bwa`;}
        print OUT "/ifs1/ST_SYSBIO/USER/liudb2/bin/TumorCare/bin/samtools index $bam;echo complete 8\n";

        print OUT "perl /ifs5/ST_COMG/USER/baiyali/RNA_editing/pan-cancer/SingleCel/human/hg19/4-parameter/Edit_detecting/bin/extract_reads_SE.pl $ba
m $output/bwa/mutation.read.fq  $output/binomFilter.filt.fre.strand.end.simplerepeat.homopolymer.dbsnp.snp.txt $output/readsvnum.txt;echo complete 9\n
";


        ###
        print OUT "$bwa aln -o 1 -e 50 -m 100000 -l 32 -k 2 -t 4 -L -i 15 $database $output/bwa/mutation.read.fq > $output/bwa/mutation.read.sai;echo 
complete 13\n";
        ###
        print OUT "$bwa samse $database $output/bwa/mutation.read.sai $output/bwa/mutation.read.fq  |/ifs1/ST_SYSBIO/USER/liudb2/bin/TumorCare/bin/sam
tools view -b -S -t $database.fai - > $output/bwa/mutation.readS.bam;echo complete 14\n";


        ###
        print OUT "perl /ifs5/ST_COMG/USER/baiyali/RNA_editing/pan-cancer/SingleCel/human/hg19/4-parameter/Edit_detecting/bin/bwa.fre.filt.SE.pl $outp
ut/bwa/mutation.readS.bam $output/binomFilter.filt.fre.strand.end.simplerepeat.homopolymer.dbsnp.snp.txt $refgene $output/readsvnum.txt $output/binomF
ilter.filt.fre.strand.end.simplerepeat.homopolymer.dbsnp.snp.bwa.txt;echo complete 15\n";
        ###
        print OUT "/ifs1/ST_SYSBIO/USER/liudb2/bin/msort -k m1 -k n2 $output/binomFilter.filt.fre.strand.end.simplerepeat.homopolymer.dbsnp.snp.bwa.tx
t > $output/binomFilter.filt.fre.strand.end.simplerepeat.homopolymer.dbsnp.snp.bwa.sort.txt;echo complete 16\n";
        ###        
        print OUT "perl /ifs1/ST_SYSBIO/USER/liudb2/RNA_editing/bladder/xiongheng/pipeline/bin/annotation.new.pl $output/binomFilter.filt.fre.strand.e
nd.simplerepeat.homopolymer.dbsnp.snp.bwa.sort.txt /ifs1/ST_SYSBIO/USER/liudb2/RNA_editing/bladder/xiongheng/pipeline/database/alu/hg19.alu.bed $outpu
t/annotation.new;echo complete 17\n";
        ###
        print OUT "perl $bin/reads.support.pl $output/annotation.new.anno.txt 2 $output/annotation.new.anno.filt.txt;echo complete 18\n";
        ###有问题
        print OUT "perl /ifs5/ST_CANCERR_COMG/USER/qiusi/pipeline/RNA_editing/strand.anno.pl $output/annotation.new.anno.filt.txt /ifs1/ST_SYSBIO/USER
/liudb2/RNA_editing/bladder/xiongheng/pipeline/bin/lincRNA/hg19_refGene.txt /ifs1/ST_SYSBIO/USER/liudb2/RNA_editing/bladder/xiongheng/pipeline/databas
e/alu/hg19.alu.bed $output/annotation.new.anno.filt.lincRNA.txt;echo complete 19\n";
       
        ###过滤父本snp
#       my $dnabam = "/ifs5/ST_COMG/USER/liwenhui/RNA_editing/ADARs/pan-cancer/pan/3.SingleCell/RNA_editing/DNA/SRR689250-1/SRR689250/brecal/SRR689250
.brecal.bam";
 #       unless (-d "$output/DNA") {`mkdir $output/DNA `;}
  #      ###有问题
#       print OUT "$bin/MismatchStat -i $dnabam -x 1000000 -o $output/DNA/dna.stat.txt -u;echo complete 20\n";
 #       #
#       print OUT "/usr/bin/perl $bin/pileupSite.pl --bam $dnabam  --reg $output/binomFilter.filt.fre.strand.end.simplerepeat.homopolymer.txt --out $o
utput/DNA/DNA.pileup --ref $chr --minMQ 20 --minBQ 20 --uniq;echo complete 21\n";
#       print OUT "$bin/pileupStat -i $output/DNA/DNA.pileup -q 20 -v $output/DNA/dna.stat.txt -o $output/DNA/dna.mut.txt.gz;echo complete 22\n";
#       print OUT "perl $bin/filterDnaMut.pl $output/DNA/dna.mut.txt.gz > $output/DNA/dna.mut.filt.txt;echo complete 23\n";
#       my $dna_mut = "$output/DNA/dna.mut.filt.txt";
#       print OUT "perl $bin/pos_filter.pl --input $output/annotation.new.anno.filt.lincRNA.txt --loci $dna_mut --out $output/annotation.new.anno.filt
.lincRNA.snp.txt --filt 2;echo complete 24\n";

    close OUT;
}
