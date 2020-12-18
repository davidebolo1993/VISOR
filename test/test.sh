#!/bin/bash

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
VISOR HACk -g smallref.fa -b hack.bed -o hack_out
echo ""
echo ""
echo "Testing VISOR SHORtS"
echo ""
echo ""
VISOR SHORtS -h
VISOR SHORtS -g smallref.fa -b shorts.bed -s hack_out -o shorts_out_singleclone --coverage 1 #this is just for testing, no need to spend a lot of time on this
VISOR SHORtS -g smallref.fa -b shorts.bed -s hack_out hack_out -o shorts_out_multiclone --coverage 1 --clonefraction 50 50 --tag #this is just for testing, no need to spend a lot of time on this
VISOR SHORtS -g smallref.fa -b shorts.bed -s hack_out -o shorts_out_strandseq --coverage 1 --strandseq --cells 1 --sce sce.bed --noise 5.0
echo ""
echo ""
echo "Testing VISOR LASeR"
echo ""
echo ""
VISOR LASeR -h
VISOR LASeR -g smallref.fa -b shorts.bed -s hack_out -o laser_out_singleclone --coverage 0.1 #this is just for testing, no need to spend a lot of time on this
VISOR LASeR -g smallref.fa -b shorts.bed -s hack_out hack_out -o laser_out_multiclone --coverage 0.1 --clonefraction 50 50 --tag #this is just for testing, no need to spend a lot of time on this
echo ""
echo ""
echo "Testing VISOR XENIA"
echo ""
echo ""
VISOR XENIA -s hack_out/ -b shorts.bed -o xenia_out --coverage 1 #this is just for testing, no need to spend a lot of time on this
echo ""
echo ""
