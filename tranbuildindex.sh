###建立hg19 index的准备文件
1、hg19.fa
2、hg19.merge.gtf 
/ifs1/ST_SINGLECELL/USER/jiangrunze/tool/tophat-2.0.12.Linux_x86_64/bin/bowtie2-build hg19.fa hg19

###建立index的main步骤
/ifs1/ST_SINGLECELL/USER/jiangrunze/tool/tophat-2.0.12.Linux_x86_64/tophat2 -G hg19.merge.gtf --transcriptome-index=bowtie2/transcriptomeindex/tranindex bowtie2/hg19/hg19

###如果是bowtie1，ti'h 替换为bowtie1
