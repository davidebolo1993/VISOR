#!/bin/bash

cd $1

find . -type f -name "*R1.fq" | sort | xargs cat > $2"."$3"_S1_L"$4"_R1_001.fastq" && find . -type f -name "*R1.fq" -delete
find . -type f -name "*R2.fq" | sort | xargs cat > $2"."$3"_S1_L"$4"_R2_001.fastq" && find . -type f -name "*R2.fq" -delete
find . -type f -name "*.fa" -delete && find . -type f -name "*.fai" -delete
