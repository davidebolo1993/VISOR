SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

echo "Creating test folder"
mkdir test && cd test

echo "Download reference"
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

echo "Download phased variants for chr22"
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi

echo "Subset to HG00732"
bcftools view -O b -o HG00732.bcf -s HG00732 -m2 -M2 -c 1 -C 1 ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
bcftools index HG00732.bcf

echo "Splitting het variants to .bed files"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' HG00732.bcf | grep "1|0" > h1.bed
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' HG00732.bcf | grep "0|1" > h2.bed

echo "Writing SNPs to VISOR .bed format"

awk 'OFS=FS="\t"''{print $1, ($2 -1), $2, "SNP", $4}' h1.bed > VISOR.h1.SNPs.bed
awk 'OFS=FS="\t"''{print $1, ($2 -1), $2, "SNP", $4}' h2.bed > VISOR.h2.SNPs.bed


echo "Subsetting reference to chr22"

samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa chr22 > chr22.fa

echo "Generating the 2 haplotypes with SNPs"

VISOR HACk -g chr22.fa -bed VISOR.h1.SNPs.bed VISOR.h2.SNPs.bed -o Templates

echo "Cloning Examples folder with .bed files for tests"

git clone https://github.com/davidebolo1993/VISOR/tree/master/Examples/test1

echo "Generating SVs in the 2 haplotypes with some random SNPs in the 2 clones"

VISOR HACk -g Templates/h1.fa -bed test1/VISOR.h1.SVs.bed -o cloneh1
VISOR HACk -g Templates/h2.fa -bed test1/VISOR.h2.SVs.bed -o cloneh2
mkdir clone
cd cloneh2 && mv h1.fa h2.fa && cd ..
mv cloneh1/h1.fa clone/
mv cloneh2/h2.fa clone/
rm -r cloneh1
rm -r cloneh2

echo "Simulating data with 80% from clone, 20% from reference"

VISOR SHORtS -g chr22.fa -s clone -bed test1/sim.bed -c 40 -o cloneout

echo "Done"
