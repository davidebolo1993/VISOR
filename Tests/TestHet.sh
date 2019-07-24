#!/bin/bash

#run in a py36 environment with VISOR installed

echo "Creating test folder"
mkdir Htest && cd Htest

echo "Downloading reference"
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

echo "Downloading phased SNPs for chr22"
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi

echo "Subsetting to HG00732"
bcftools view -O b -o HG00732.bcf -s HG00732 -m2 -M2 -c 1 -C 1 ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
bcftools index HG00732.bcf

echo "Splitting het variants into different BED"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' HG00732.bcf | grep "1|0" > h1.bed
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' HG00732.bcf | grep "0|1" > h2.bed
cat h1.bed h2.bed | awk '{print $2}' | sort > allSNPs.txt

echo "Writing SNPs to BED for HACk"

awk 'OFS=FS="\t"''{print $1, ($2 -1), $2, "SNP", $4, "0"}' h1.bed > VISOR.h1.SNPs.bed
awk 'OFS=FS="\t"''{print $1, ($2 -1), $2, "SNP", $4, "0"}' h2.bed > VISOR.h2.SNPs.bed

echo "Subsetting reference to chr22"

samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa chr22 > chr22.fa

echo "Generating the 2 FASTA haplotypes with SNPs"

VISOR HACk -g chr22.fa -bed VISOR.h1.SNPs.bed VISOR.h2.SNPs.bed -o Templates

echo "Generating HACk.bed with variants"

mkdir files && cd files

echo -e "chr22\t20000000\t22000000\tdeletion\tNone\t0\nchr22\t35000000\t37000000\ttandem duplication\t2\t0" > VISOR.clone1.h1.SVs.bed
echo -e "chr22\t20000000\t22000000\tdeletion\tNone\t0\nchr22\t42000000\t45000000\tdeletion\tNone\t0" > VISOR.clone2.h1.SVs.bed
echo -e "chr22\t0\t50818468\t100.0\t100.0" > VISOR.sim.bed

cd ..

wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestHet/pileup2base.pl
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestHet/plotBAFCOV.R

echo "Inserting SVs in the 2 clones"

VISOR HACk -g Templates/h1.fa -bed files/VISOR.clone1.h1.SVs.bed -o clone1
cp chr22.fa clone1/ && mv clone1/chr22.fa clone1/h2.fa

VISOR HACk -g Templates/h1.fa -bed files/VISOR.clone2.h1.SVs.bed -o clone2
cp chr22.fa clone2/ && mv clone2/chr22.fa clone2/h2.fa

echo "Simulating data. Clone 1: 65%; Clone 2: 30 %; Reference: 5%"

VISOR SHORtS -g chr22.fa -s clone1/ clone2/ Templates/ -bed files/VISOR.sim.bed -c 150 -o cloneout --clonefraction 65.0 30.0 5.0 -- threads 7

echo "Running mpileup"

samtools mpileup -f GRCh38_full_analysis_set_plus_decoy_hla.fa -r chr22 cloneout/sim.srt.bam > cloneout/tumor.mpileup
grep -f allSNPs.txt cloneout/tumor.mpileup > cloneout/het.tumor.mpileup

echo "Extracting info"

perl pileup2base.pl cloneout/het.tumor.mpileup 20 cloneout/het.tumor.pileup2base

echo "Plotting"

R --slave --args cloneout/het.tumor.pileup2base,cloneout/testhet.pdf < plotBAFCOV.R

echo "PDF saved in cloneout/testhet.pdf"

echo "Done"
