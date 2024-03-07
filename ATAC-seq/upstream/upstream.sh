#!/bin/bash
for i in *.bam
do 
macs2 callpeak -f BAM -t $i -n $i --outdir ./macs2/ --bdg -q 0.05
done