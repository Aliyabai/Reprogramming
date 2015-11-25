#!usr/bin.perl -w
use strict;
die "Usage: perl $0 <SRR.list>\n" if (@ARGV != 1);
my $maindir = "/ifs5/ST_COMG/USER/baiyali/RNA_editing/Reprogramming/data_process";
my $work = "/ifs5/ST_COMG/USER/baiyali/RNA_editing/Reprogramming/data_process/work/bamprocess";
my @samplefiles = `ls $maindir/SRR*/tophat/accepted_hits.bam`;
chomp @samplefiles;
my $markDuplicatesdir = "java -Xmx4g -jar /ifs1/ST_SYSBIO/USER/liudb2/bin/picard/picard-tools-1.84/MarkDuplicates.jar";
my $reorder = "java -Xmx4g -jar /ifs1/ST_SYSBIO/USER/liudb2/bin/picard/picard-tools-1.84/ReorderSam.jar";
my $rmdupparameters = " REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true";
my $reorderparameters = "REFERENCE=/ifs1/ST_SYSBIO/USER/liudb2/RNA_editing/bladder/xiongheng/database/hg19index/hg19.fa VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true";
my $baseRecalibratorparameters = "-nct 3 -U ALLOW_N_CIGAR_READS";
my $ref = "-R /ifs1/ST_SYSBIO/USER/liudb2/RNA_editing/bladder/xiongheng/database/hg19index/hg19.fa";
my $databasedir = "/ifs1/ST_SYSBIO/USER/liudb2/database/GATK/2.8/hg19";
my $javadir1 = "/ifs1/ST_SYSBIO/USER/liudb2/bin/java/jre1.7.0/bin/java -Xmx6g -jar /ifs1/ST_SYSBIO/USER/liudb2/bin/GATK/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar -T BaseRecalibrator -l INFO";
my $javadir2 = "/ifs1/ST_SYSBIO/USER/liudb2/bin/java/jre1.7.0/bin/java -Xmx6g -jar /ifs1/ST_SYSBIO/USER/liudb2/bin/GATK/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar -T PrintReads -l INFO";

open IN, "<$ARGV[0]" or die "$!";  
my @indexline=<IN>;
chomp @indexline;
close IN;
#print "@indexline\n";

for my $indexline (@indexline){
    chomp $indexline;# print "$indexline\n";
    for my $eachfile (@samplefiles){
        if ($eachfile =~/$indexline/){
           unless (-d "$maindir/$indexline/bamprocess/"){`mkdir $maindir/$indexline/bamprocess/ `;}
           #print "id:$indexline\tfile:$eachfile\n";
           open OUT, ">$work/bam.process.$indexline.sh" or die $!;
           my $sample = $indexline; 
           print OUT "$reorder INPUT=$eachfile OUTPUT=$maindir/$sample/bamprocess/$sample.reorder.bam $reorderparameters;echo 1 complete\n";
           print OUT "$markDuplicatesdir INPUT=$maindir/$sample/bamprocess/$sample.reorder.bam OUTPUT=$maindir/$sample/bamprocess/$sample.reorder.mark
dup.bam METRICS_FILE=$maindir/$sample/bamprocess/$sample.bam.rmdup.met $rmdupparameters;echo 2 complete\n";
           print OUT "$javadir1 -I $maindir/$sample/bamprocess/$sample.reorder.markdup.bam $ref -knownSites $databasedir/dbsnp_138.hg19.vcf.gz -knownS
ites $databasedir/1000G_phase1.indels.hg19.vcf.gz -knownSites $databasedir/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz -o $maindir/$sample/bampro
cess/$sample.reorder.markdup.baseRecal.grp $baseRecalibratorparameters;echo 3 complete\n";
           print OUT "$javadir2 -I $maindir/$sample/bamprocess/$sample.reorder.markdup.bam $ref -BQSR $maindir/$sample/bamprocess/$sample.reorder.mark
dup.baseRecal.grp -U ALLOW_N_CIGAR_READS -o $maindir/$sample/bamprocess/$sample.reorder.markdup.baseRecal.bam;echo 4 complete\n";
           close OUT;
        }        
    }
}
