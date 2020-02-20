#!/bin/bash

cd $1
awk -v var="$2" '{if(NR%4==1) $0=sprintf("@"var"_%d",(1+i++)); if (NR%4==3) $0=sprintf("+"var"_%d",(1+l++)); print;}' sim_0001.fastq > tmp.fq
mv tmp.fq sim_0001.fastq