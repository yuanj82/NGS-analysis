#!/bin/bash

set -e
set -o pipefail

############# environmental preparation #############

conda activate RNA-seq # 使用conda环境

mkdir project
wkd='./project/'
projectName='project'
cd $wkd
mkdir $projectName && cd $projectName

# SRR_Acc_List.txt文件放在project根目录便于调用
mkdir sra fastq_gz fastqc_reports clean sorted genome_index compared

############# data preparation #############

nohup prefetch -O . $(<SRR_Acc_List.txt) &
wait # 等待prefetch完成

cat ./SRR_Acc_List.txt | while read id; do
    mv -f "${id}/${id}.sra" ./sra
    rm -rf "${id}"
done

############# sra to fastq #############

ls *.sra | parallel fastq-dump --gzip --split-3 {} --outdir ./fastqgz

############# FastQc quality control #############

ls ./fastqgz/*.fastq.gz | xargs fastqc -t 2 -O ./fastqc_reports
multiqc ./fastqc_report  # 合并质量检测报告

############# Trim_galore filtering sequence #############

cat ./SRR_Acc_List.txt | while read id; do 
    nohup trim_galore \
    -q 20 \
    --length 36 \
    --max_n 3 \
    --stringency 3 \
    --fastqc \
    --paired \
    -o ./clean/ \
    ./fastqgz/${id}_1.fastq.gz ./fastqgz/${id}_2.fastq.gz &
done
wait # 等待trim_galore完成

############# Reply to reference genome #############

genome='../oryza_sativa.fa'   # 基因组文件路径（需提前下载解压）
species='oryza_sativa'  # 物种名称（作为索引名称）
cd ./genome_index
cp ${genome} ./

hisat2-build -p 4 ${species}.fa ${species}
cd ..

cat ./SRR_Acc_List.txt | while read id;
do
    hisat2 \
    -x ./genome_index/${species} \
    -p 5 \
    -1 ./clean/${id}_1.fastq.gz \
    -2 ./clean/${id}_2.fastq.gz \
    -S ./compared/${id}.sam
done

############# sam to bam #############

cat ./SRR_Acc_List.txt | while read id ;do
    samtools \
    sort \
    -n -@ 5 \
    ./compared/${id}.sam \
    -o ./sorted/${id}.bam
done

############# featureCounts #####

cd ./sorted
bam=$(ls *)
gtf='../oryza_sativa.gff3'
featureCounts \
-T 5 \
-t exon \
-g Name \
-a ${gtf} \
-o gene_name.counts \
-p ${bam}