gzip -dc /ifs1/pub/database/ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz > /ifs5/ST_COMG/USER/baiyali/database/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa;echo 1 complete
gzip -dc /ifs1/pub/database/ftp.ensembl.org/pub/release-81/gtf/homo_sapiens/Homo_sapiens.GRCh38.81.gtf.gz > /ifs5/ST_COMG/USER/baiyali/database/hg38/Homo_sapiens.GRCh38.81.gtf;echo 2 complete
ln -s Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa hg38.fa;echo 3 complete
/home/xiongheng/bin/bowtie2-2.2.3/bowtie2-build hg38.fa hg38;echo 4 complete
samtools faidx *.fa
