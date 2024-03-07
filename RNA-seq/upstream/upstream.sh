#!/bin/bash

############# environmental preparation #############

conda env create -f env.yml # 建立conda环境
conda activate genomics # 使用conda环境

mkdir project
wkd='./project/'
projectName='project'
cd $wkd
mkdir $projectName && cd $projectName
 
mkdir sra fastq_gz fastqc_reports clean sorted 

############# data preparation #############

cd sra

nohup prefetch -O . $(<SRR_Acc_List.txt) &
# cat SRR_Acc_List.txt | while read id; do (prefetch ${id});done
ps -ef | grep prefetch | awk '{print $2}' | while read id; do kill ${id}; done ## kill download processes
 
cat SRR_Acc_List.txt | while read id; do mv -f $id/$id.sra  ./; done 
cat SRR_Acc_List.txt | while read id; do rm -rf $id; done #分两步走，检出全部移动完毕再删除文件夹
 
############# sra to fastq #############

for i in *sra
do
	fastq-dump --gzip --split-3 $i -f ../fastq_gz
done

############# FastQc quality control #############

cd ../fastq_gz
ls *.fastq.gz | xargs fastqc -t 2 -O ../fastqc_reports
cd ../fastqc_reports && multiqc ./
    
############# Trimmatic filtering sequence #############

cd ../fastq_gz

fastq_1=($(find . -type f -name "*_1.fastq.gz"))
fastq_2=($(find . -type f -name "*_2.fastq.gz"))

for ((i=0; i<${#fastq_1[@]}; i++)); do
    name=$(basename ${fastq_1[$i]} _1.fastq.gz)
    trimmomatic PE -threads 4 -phred33 "${fastq_1[$i]}" "${fastq_2[$i]}" -summary "${name}.summary" -baseout "${name}.fastq.gz" LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 HEADCROP:15 MINLEN:36
done

############# sam to bam #############

hisat2-build -p 4 oryza_sativa.fa oryza_sativa
 
############# sam to bam #############
cd ../aligned
ls *.sam| while read id ;do (samtools sort -O bam -@ 5 -o $(basename ${id} ".sam").bam ${id});done
rm *.sam
# 为bam文件建立索引
ls *.bam |xargs -i samtools index {}
# 比对结果统计
ls *.bam |while read id ;do ( samtools flagstat -@ 1 $id >  $(basename ${id} ".bam").flagstat  );done
 
############# featureCounts #####
gtf='/home/fleame/public/references/gtf/gencode.v32.annotation.gtf.gz'   
featureCounts -T 5 -p -t exon -g gene_id  -a $gtf -o  all.id.txt  *.bam  1>counts.id.log 2>&1 &
