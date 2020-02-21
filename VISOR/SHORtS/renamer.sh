#!/bin/bash

cd $1
awk -v var="$2" '{if(NR%4==1) $0=sprintf("@"var"_%d_"substr($0,2),(1+i++)); print;}' region.1.fq > region.1.tmp.fq
awk -v var="$2" '{if(NR%4==1) $0=sprintf("@"var"_%d_"substr($0,2),(1+i++)); print;}' region.2.fq > region.2.tmp.fq
mv region.1.tmp.fq region.1.fq
mv region.2.tmp.fq region.2.fq