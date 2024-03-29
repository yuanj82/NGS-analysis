#!/bin/bash

set -e

############# environmental preparation #############

conda activate RNA-seq # 使用conda环境

mkdir project
wkd='./project/'
projectName='project'
cd $wkd
mkdir $projectName && cd $projectName

mkdir sra fastq_gz fastqc_reports clean sorted genome_index compared

############# data preparation #############

cd sra

nohup prefetch -O . $(<SRR_Acc_List.txt) &  # 将进程挂在后台进行

cat SRR_Acc_List.txt | while read id; do
    mv -f "${id}/${id}.sra" ./
    rm -rf "${id}"
done
 
############# sra to fastq #############

ls *.sra | parallel fastq-dump --gzip --split-3 {} --outdir ../fastqgz

############# FastQc quality control #############

cd ../fastq_gz
ls *.fastq.gz | xargs fastqc -t 2 -O ../fastqc_reports
cd ../fastqc_reports && multiqc ./  # 合并质量检测报告
    
############# Trim_galore filtering sequence #############

cd ../fastq_gz

cat SRR_Acc_List.txt | while read id; do 
    nohup trim_galore \
    -q 20 \
    --length 36 \
    --max_n 3 \
    --stringency 3 \
    --fastqc \
    --paired \
    -o ../clean/ \
    ${id}_1.fastq.gz ${id}_2.fastq.gz & 
done

############# Reply to reference genome #############

genome='../oryza_sativa.fa'
species='oryza_sativa'
cd ../genome_index

nohup hisat2-build -p 4 ${genome} ${species} &

less ../SRR_Acc_List.txt | while read id;
do
    hisat2 \
    -x ${species}/${species} \
    -p 5 \
    -1 ${id}_1.fastq.gz \
    -2 ${id}_2.fastq.gz \
    -S ${id}.sam
done

mv *.sam ../compared

############# sam to bam #############

cd ../compared
ls *.sam | while read id ;do
    samtools \
    sort \
    -n -@ 5 \
    ${id}.sam \
    -o ${id}.bam
done

rm *.sam
mv ${id}.bam ../sorted
 
############# featureCounts #####

cd ../sorted
bam=$(ls *)
gtf='../oryza_sativa.gff3'
featureCounts \
-T 5 \
-t exon \
-g Name \
-a ${gtf} \
-o gene_name.counts \
-p ${bam}


