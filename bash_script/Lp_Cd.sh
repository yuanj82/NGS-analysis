############# Environmental preparation #############

conda create -n rna-seq python=3.7

conda activate rna-seq

conda install fastqc trimmomatic hisat2 samtools subread #subread 即 featureCounts

wkd='./project/' # 指定工作目录
projectName='Lp-Cd' # 指定项目名称
cd $wkd
mkdir $projectName && cd $projectName # 创建项目文件夹
 
mkdir fastq_gz fastqc_reports clean sorted hisat2 # 创建所需的文件夹

############# Data download #############

cd fastq_gz

# 下载测序数据
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C0_1.R1.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C0_1.R2.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C0_2.R1.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C0_2.R2.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C0_3.R1.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C0_3.R2.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C500_1.R1.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C500_1.R2.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C500_2.R1.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C500_2.R2.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C500_3.R1.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C500_3.R2.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C50_1.R1.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C50_1.R2.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C50_2.R1.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C50_2.R2.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C50_3.R1.fastq.gz
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lp/C50_3.R2.fastq.gz

rm wget-log # 移除下载日志

############# Quality control #############

fastqc *.fastq.gz -o ../fastqc_reports # 进行质量检测，将报告存放到fastqc_reports文件夹
cd ../fastqc_reports && multiqc ./ # 合并检测报告

cd ../fastq_gz

############# Sequence filtering #############

# 进行序列过滤
for i in {0,50,500}
do
    for j in {1,2,3}
    do
        trimmomatic PE -threads 4 -phred33 C${i}_${j}.R1.fastq.gz C${i}_${j}.R2.fastq.gz -summary Lp_${i}.summary -baseout C${i}_${j}.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 HEADCROP:15 MINLEN:36
    done
done

# 将过滤后的数据移动到clean文件夹
mv *U.fastq.gz ../clean
mv *P.fastq.gz ../clean

############# Hisat2 comparison #############

# 下载黑麦草参考基因组
wget https://data-1310575002.cos.ap-chengdu.myqcloud.com/Lolium_perenne.MPB_Lper_Kyuss_1697.dna.toplevel.fa

rm wget-log

# 重命名基因组序列文件
mv Lolium_perenne.MPB_Lper_Kyuss_1697.dna.toplevel.fa Lolium_perenne.fa

mkdir index # 创建文件夹用于存放hisat2索引

hisat2-build -p 4 Lolium_perenne.fa Lolium_perenne # 建立索引

# 将索引移动到index文件夹
mv *.ht2 index/
mv Lolium_perenne.fa index/

# 将reads比对到参考基因组
for i in {0,50,500}
do
    for j in {1,2,3}
    do
        hisat2 -x index/Lolium_perenne -p 5 -1 C${i}_${j}_1P.fastq.gz -2 C${i}_${j}_2P.fastq.gz -S C${i}_${j}.sam
    done
done

# 将hisat2相关的文件归档到hisat2文件夹
mv index ../hisat2
mv *.sam ../hisat2

############# Sort compression #############

cd ../hisat2

# 进行排序压缩
for i in {0,50,500}
do
    for j in {1,2,3}
    do
        samtools sort -n -@ 5 C${i}_${j}.sam -o C${i}_${j}.bam
    done
done

# 将排序压缩后的文件归档到sorted文件夹
mv *.bam ../sorted

############# Counting statistics #############

cd ../sorted

# 下载黑麦草注释文件
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/gff3/lolium_perenne/Lolium_perenne.MPB_Lper_Kyuss_1697.57.gff3.gz

rm wget-log

# 解压注释文件
gzip -d Lolium_perenne.MPB_Lper_Kyuss_1697.57.gff3.gz

mv Lolium_perenne.MPB_Lper_Kyuss_1697.57.gff3 Lolium_perenne.gff3

# 创建数组，存储所有bam文件
bam_files=(*.bam)

# 进行计数统计
if [ ${#bam_files[@]} -gt 0 ]; then
    featureCounts -T 5 -t exon -g Name -a Lolium_perenne.gff3 -o  gene.counts -p "${bam_files[@]}"
else
    echo "No BAM files found in the current directory."
fi

mv gene.counts* ../