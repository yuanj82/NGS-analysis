
nohup prefetch -O . $(<SRR_Acc_List.txt) &

nohup ls *.fastq.gz | xargs fastqc -t 2 -O ../fastqc_reports &


for i in {57..62}
do
    trim_galore -q 20 --length 36 --max_n 3 --stringency 3 --fastqc --paired -o SRR58746{i}_1_trim.fastq.gz SRR58746{i}_2_trim.fastq.gz SRR58746{i}_1.fastq.gz SRR58746{i}_2.fastq.gz > {log} 2>&1 
done



nohup wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-51/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz &

#!/bin/bash
for i in *.bam
do 
macs2 callpeak -f BAM -t $i -n $i --outdir ./macs2/ --bdg -q 0.05
done