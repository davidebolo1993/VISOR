#!/bin/bash

cd $1

grep -f watsonreads.txt -A 3 region.1.fq | grep -v '^--' > watson.1.fq
grep -f watsonreads.txt -A 3 region.2.fq | grep -v '^--' > watson.2.fq
grep -f crickreads.txt -A 3 region.1.fq | grep -v '^--' > crick.1.fq
grep -f crickreads.txt -A 3 region.2.fq | grep -v '^--' > crick.2.fq
rm region.1.fq && rm region.2.fq
