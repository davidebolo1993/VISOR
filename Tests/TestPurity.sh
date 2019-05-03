#!/bin/bash

#run in a py36 environment with VISOR installed

echo "Create test folder"
mkdir Ptest && cd Ptest

echo "Download reference"
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

echo "Download phased variants for chr22"
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi

echo "Subset to HG00732"
bcftools view -O b -o HG00732.bcf -s HG00732 -m2 -M2 -c 1 -C 1 ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
bcftools index HG00732.bcf

echo "Split het variants to different BED"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' HG00732.bcf | grep "1|0" > h1.bed
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' HG00732.bcf | grep "0|1" > h2.bed
cat h1.bed h2.bed | awk '{print $2}' | sort > allSNPs.txt

echo "Write SNPs to BED for HACk"

awk 'OFS=FS="\t"''{print $1, ($2 -1), $2, "SNP", $4}' h1.bed > VISOR.h1.SNPs.bed
awk 'OFS=FS="\t"''{print $1, ($2 -1), $2, "SNP", $4}' h2.bed > VISOR.h2.SNPs.bed

echo "Subset reference to chr22"

samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa chr22 > chr22.fa

echo "Generate the 2 FASTA haplotypes with SNPs"

VISOR HACk -g chr22.fa -bed VISOR.h1.SNPs.bed VISOR.h2.SNPs.bed -o Templates

echo "Get example folder with BED files for tests from VISOR git"

mkdir files && cd files

wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestPurity/VISOR.h1.SVs.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestPurity/VISOR.h2.SVs.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestPurity/VISOR.sim.bed

cd ..

wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestPurity/pileup2base.pl
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestPurity/plotBAFCOV.R

echo "Generate SVs in clone with some random SNPs"

VISOR HACk -g Templates/h1.fa -bed files/VISOR.h1.SVs.bed -o cloneh1
VISOR HACk -g Templates/h2.fa -bed files/VISOR.h2.SVs.bed -o cloneh2
mkdir clone
cd cloneh2 && mv h1.fa h2.fa && cd ..
mv cloneh1/h1.fa clone/
mv cloneh2/h2.fa clone/
rm -r cloneh1
rm -r cloneh2

echo "Simulate data. Clone: 80%; Reference contamination: 20%"

VISOR SHORtS -g chr22.fa -s clone/ Templates/ -bed files/VISOR.sim.bed -c 100 -o cloneout -cf 80.0 20.0

echo "Simulations done"

echo "Run mpileup"

samtools mpileup -f GRCh38_full_analysis_set_plus_decoy_hla.fa -r chr22 cloneout/sim.srt.bam > cloneout/tumor.mpileup
grep -f allSNPs.txt cloneout/tumor.mpileup > cloneout/het.tumor.mpileup

echo "Done"

echo "Extract info"

perl pileup2base.pl cloneout/het.tumor.mpileup 20 cloneout/het.tumor.pileup2base

echo "Done"

echo "Plot"

R --slave --args cloneout/het.tumor.pileup2base,cloneout/testpurity.pdf < plotBAFCOV.R

echo "PDF saved in cloneout/testhet.pdf"
echo "Done"
