#!/bin/bash
echo "Checking executable installation"
echo ""
echo ""
samtools --version
minimap2 --version
bedtools --version
badread --version
echo ""
echo ""
echo "Checking VISOR installation"
echo ""
echo ""
VISOR -h
echo ""
echo ""
echo "Testing VISOR HACk"
echo ""
echo ""
while true;do
wget -T 15 -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa && break
done
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa chr1 > smallref.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa chr2 >> smallref.fa
VISOR HACk -h
VISOR HACk -g smallref.fa -b test/hack.bed -o test/hack_out
echo ""
echo ""
echo "Testing VISOR SHORtS"
echo ""
echo ""
VISOR SHORtS -h
VISOR SHORtS -g smallref.fa -b test/shorts.bed -s test/hack_out -o test/shorts_out_singleclone --coverage 0.0001 #this is just for testing, no need to spend a lot of time on this
VISOR SHORtS -g smallref.fa -b test/shorts.bed -s test/hack_out test/hack_out -o test/shorts_out_multiclone --coverage 0.0001 --clonefraction 50 50 --tag #this is just for testing, no need to spend a lot of time on this
VISOR SHORtS -g smallref.fa -b test/shorts.bed -s test/hack_out -o test/shorts_out_strandseq --coverage 0.0001 --strandseq --cells 1 --sce test/sce.bed --noise 5.0
echo ""
echo ""
echo "Testing VISOR LASeR"
echo ""
echo ""
VISOR LASeR -h
VISOR LASeR -g smallref.fa -b test/shorts.bed -s test/hack_out -o test/laser_out_singleclone --coverage 0.00001 --threads 2 #this is just for testing, no need to spend a lot of time on this
VISOR LASeR -g smallref.fa -b test/shorts.bed -s test/hack_out test/hack_out -o test/laser_out_multiclone --coverage 0.00001 --clonefraction 50 50 --tag --threads 2 #this is just for testing, no need to spend a lot of time on this
echo ""
echo ""
echo "Testing VISOR XENIA"
echo ""
echo ""
VISOR XENIA -s test/hack_out -b test/shorts.bed -o test/xenia_out --coverage 0.0001 #this is just for testing, no need to spend a lot of time on this
echo ""
echo ""
