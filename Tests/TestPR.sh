#!/bin/bash

#run in a py36 environment with VISOR installed

echo "Creating test folder"
mkdir PRtest && cd PRtest

echo "Downloading reference"

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

echo "Downloading phased SVs"

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/GRCh38_positions/ALL.wgs.integrated_sv_map_v2_GRCh38.20130502.svs.genotypes.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/GRCh38_positions/ALL.wgs.integrated_sv_map_v2_GRCh38.20130502.svs.genotypes.vcf.gz.tbi


echo "Subsetting to HG00732"


bcftools view -O b -o HG00732.bcf -s HG00732 -m2 -M2 -c 1 -C 1 ALL.wgs.integrated_sv_map_v2_GRCh38.20130502.svs.genotypes.vcf.gz
bcftools index HG00732.bcf



echo "Getting DELETIONS, INSERTIONS and INVERSIONS and write .bed for VISOR"

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE[\t%SAMPLE=%GT]\n' HG00732.bcf  | grep -w "DEL"  | grep "1|0" | awk 'OFS=FS="\t"''{if ($3 -$2 >1000) print $1, $2, $3, "deletion", "None", "0"}' | sed 's/^/chr/' > VISOR.SV.h1.HG00732.bed
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE[\t%SAMPLE=%GT]\n' HG00732.bcf  | grep -w "DEL"  | grep "0|1" | awk 'OFS=FS="\t"''{if ($3 -$2 >1000) print $1, $2, $3, "deletion", "None", "0"}' | sed 's/^/chr/' > VISOR.SV.h2.HG00732.bed

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE[\t%SAMPLE=%GT]\n' HG00732.bcf  | grep -w "DUP" | grep "1|0" | awk 'OFS=FS="\t"''{if ($3 -$2 >1000) print $1, $2, $3, "tandem duplication", "2", "0"}' | sed 's/^/chr/' >> VISOR.SV.h1.HG00732.bed
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE[\t%SAMPLE=%GT]\n' HG00732.bcf  | grep -w "DUP" | grep "0|1" | awk 'OFS=FS="\t"''{if ($3 -$2 >1000) print $1, $2, $3, "tandem duplication", "2", "0"}' | sed 's/^/chr/' >> VISOR.SV.h2.HG00732.bed


bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE[\t%SAMPLE=%GT]\n' HG00732.bcf  | grep -w "INV" | grep "1|0" | awk 'OFS=FS="\t"''{if ($3 -$2 >1000) print $1, $2, $3, "inversion", "None", "0"}' | sed 's/^/chr/' >> VISOR.SV.h1.HG00732.bed
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE[\t%SAMPLE=%GT]\n' HG00732.bcf  | grep -w "INV" | grep "0|1" | awk 'OFS=FS="\t"''{if ($3 -$2 >1000) print $1, $2, $3, "inversion", "None", "0"}' | sed 's/^/chr/' >> VISOR.SV.h1.HG00732.bed


echo "Generating haplotypes with SVs"

VISOR HACk -g GRCh38_full_analysis_set_plus_decoy_hla.fa -bed VISOR.SV.h1.HG00732.bed VISOR.SV.h2.HG00732.bed -o Templates


echo "Getting .bed for simulations with correct chormosomes sizes from VISOR git"

mkdir files && cd files

wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestPR/VISOR.sim.bed

cd ..


echo "Simulating short reads with different coverages"

VISOR SHORtS -g GRCh38_full_analysis_set_plus_decoy_hla.fa -s Templates/ -bed files/VISOR.sim.bed -c 10 -o cov10s -th 7 --noaddtag
VISOR SHORtS -g GRCh38_full_analysis_set_plus_decoy_hla.fa -s Templates/ -bed files/VISOR.sim.bed -c 20 -o cov20s -th 7 --noaddtag
VISOR SHORtS -g GRCh38_full_analysis_set_plus_decoy_hla.fa -s Templates/ -bed files/VISOR.sim.bed -c 30 -o cov30s -th 7 --noaddtag



echo "Running DELLY: will fail if DELLY is not installed"


delly call -g GRCh38_full_analysis_set_plus_decoy_hla.fa cov10s/sim.srt.bam -o cov10s/delly.vcf
delly call -g GRCh38_full_analysis_set_plus_decoy_hla.fa cov20s/sim.srt.bam -o cov20s/delly.vcf
delly call -g GRCh38_full_analysis_set_plus_decoy_hla.fa cov30s/sim.srt.bam -o cov30s/delly.vcf



echo "Simulating long reads with different coverages"

VISOR LASeR -g GRCh38_full_analysis_set_plus_decoy_hla.fa -s Templates/ -bed files/VISOR.sim.bed -c 10 -o cov10l -th 7 --noaddtag
VISOR LASeR -g GRCh38_full_analysis_set_plus_decoy_hla.fa -s Templates/ -bed files/VISOR.sim.bed -c 20 -o cov20l -th 7 --noaddtag
VISOR LASeR -g GRCh38_full_analysis_set_plus_decoy_hla.fa -s Templates/ -bed files/VISOR.sim.bed -c 30 -o cov30l -th 7 --noaddtag



echo "Running SVIM: will fail if SVIM is not installed"

svim alignment cov10l cov10l/sim.srt.bam
svim alignment cov20l cov20l/sim.srt.bam
svim alignment cov30l cov30l/sim.srt.bam



cat VISOR.SV.h1.HG00732.bed VISOR.SV.h2.HG00732.bed  | sortBed  | mergeBed > ALL.SVs.bed #merge overlapping SVs


echo "Generating true positives (TP), false positives (FP) and false negatives (FN) .bed files to make you easier to compute precision and recall :) "


cd cov10s

bcftools view -f "PASS" -e IMPRECISE delly.vcf | bcftools view -i MAPQ=60 -O b > delly.fltrd.bcf
bcftools index delly.fltrd.bcf
bcftools query -f '%CHROM\t%POS\t%INFO/END\n' delly.fltrd.bcf > variants.bed

bedtools intersect -a variants.bed -b ../ALL.SVs.bed -wa > TP.bed #VARIANTS CALLED BY DELLY AND MATCHING THE HIGH-CONFIDENCE SET
bedtools intersect -a variants.bed -b ../ALL.SVs.bed -v > FP.bed #VARIANTS CALLED BY DELLY AND NOT MATCHING THE HIGH-CONFIDENCE SET
bedtools intersect -a ../ALL.SVs.bed -b variants.bed -v > FN.bed #VARIANTS IN HIGH-CONFIDENCE SET NOT CALLED BY DELLY

cd ..


cd cov10l


awk '{ if($1 ~ /^#/) { print $0 } else { if($6>=3) { print $0 } } }' final_results.vcf  > svim.fltrd.vcf #filtered using the script here: https://github.com/eldariont/svim/wiki, changing threshold for low coverage data (COV < 30)
bcftools view -h svim.fltrd.vcf > header.vcf
bcftools view -H svim.fltrd.vcf | grep -v "INS" >> header.vcf
mv header.vcf svim.fltrd.vcf && bcftools view -O b svim.fltrd.vcf > svim.fltrd.bcf
bcftools index svim.fltrd.bcf
bcftools query -f '%CHROM\t%POS\t%INFO/END\n' svim.fltrd.bcf > variants.bed

bedtools intersect -a variants.bed -b ../ALL.SVs.bed -wa > TP.bed #VARIANTS CALLED BY SVIM AND MATCHING THE HIGH-CONFIDENCE SET
bedtools intersect -a variants.bed -b ../ALL.SVs.bed -v > FP.bed #VARIANTS CALLED BY SVIM AND NOT MATCHING THE HIGH-CONFIDENCE SET
bedtools intersect -a ../ALL.SVs.bed -b variants.bed -v > FN.bed #VARIANTS IN HIGH-CONFIDENCE SET NOT CALLED BY SVIM


cd ..


cd cov20s

bcftools view -f "PASS" -e IMPRECISE delly.vcf | bcftools view -i MAPQ=60 -O b > delly.fltrd.bcf
bcftools index delly.fltrd.bcf
bcftools query -f '%CHROM\t%POS\t%INFO/END\n' delly.fltrd.bcf > variants.bed

bedtools intersect -a variants.bed -b ../ALL.SVs.bed -wa > TP.bed #VARIANTS CALLED BY DELLY AND MATCHING THE HIGH-CONFIDENCE SET
bedtools intersect -a variants.bed -b ../ALL.SVs.bed -v > FP.bed #VARIANTS CALLED BY DELLY AND NOT MATCHING THE HIGH-CONFIDENCE SET
bedtools intersect -a ../ALL.SVs.bed -b variants.bed -v > FN.bed #VARIANTS IN HIGH-CONFIDENCE SET NOT CALLED BY DELLY

cd ..


cd cov20l


awk '{ if($1 ~ /^#/) { print $0 } else { if($6>=6) { print $0 } } }' final_results.vcf  > svim.fltrd.vcf #filtered using the script here: https://github.com/eldariont/svim/wiki, changing threshold for low coverage data (COV < 30)
bcftools view -h svim.fltrd.vcf > header.vcf
bcftools view -H svim.fltrd.vcf | grep -v "INS" >> header.vcf
mv header.vcf svim.fltrd.vcf && bcftools view -O b svim.fltrd.vcf > svim.fltrd.bcf
bcftools index svim.fltrd.bcf
bcftools query -f '%CHROM\t%POS\t%INFO/END\n' svim.fltrd.bcf > variants.bed

bedtools intersect -a variants.bed -b ../ALL.SVs.bed -wa > TP.bed #VARIANTS CALLED BY SVIM AND MATCHING THE HIGH-CONFIDENCE SET
bedtools intersect -a variants.bed -b ../ALL.SVs.bed -v > FP.bed #VARIANTS CALLED BY SVIM AND NOT MATCHING THE HIGH-CONFIDENCE SET
bedtools intersect -a ../ALL.SVs.bed -b variants.bed -v > FN.bed #VARIANTS IN HIGH-CONFIDENCE SET NOT CALLED BY SVIM


cd ..


cd cov30s

bcftools view -f "PASS" -e IMPRECISE delly.vcf | bcftools view -i MAPQ=60 -O b > delly.fltrd.bcf
bcftools index delly.fltrd.bcf
bcftools query -f '%CHROM\t%POS\t%INFO/END\n' delly.fltrd.bcf > variants.bed

bedtools intersect -a variants.bed -b ../ALL.SVs.bed -wa > TP.bed #VARIANTS CALLED BY DELLY AND MATCHING THE HIGH-CONFIDENCE SET
bedtools intersect -a variants.bed -b ../ALL.SVs.bed -v > FP.bed #VARIANTS CALLED BY DELLY AND NOT MATCHING THE HIGH-CONFIDENCE SET
bedtools intersect -a ../ALL.SVs.bed -b variants.bed -v > FN.bed #VARIANTS IN HIGH-CONFIDENCE SET NOT CALLED BY DELLY

cd ..


cd cov30l


awk '{ if($1 ~ /^#/) { print $0 } else { if($6>=9) { print $0 } } }' final_results.vcf  > svim.fltrd.vcf #filtered using the script here: https://github.com/eldariont/svim/wiki, changing threshold for low coverage data (COV < 30)
bcftools view -h svim.fltrd.vcf > header.vcf
bcftools view -H svim.fltrd.vcf | grep -v "INS" >> header.vcf
mv header.vcf svim.fltrd.vcf && bcftools view -O b svim.fltrd.vcf > svim.fltrd.bcf
bcftools index svim.fltrd.bcf
bcftools query -f '%CHROM\t%POS\t%INFO/END\n' svim.fltrd.bcf > variants.bed

bedtools intersect -a variants.bed -b ../ALL.SVs.bed -wa > TP.bed #VARIANTS CALLED BY SVIM AND MATCHING THE HIGH-CONFIDENCE SET
bedtools intersect -a variants.bed -b ../ALL.SVs.bed -v > FP.bed #VARIANTS CALLED BY SVIM AND NOT MATCHING THE HIGH-CONFIDENCE SET
bedtools intersect -a ../ALL.SVs.bed -b variants.bed -v > FN.bed #VARIANTS IN HIGH-CONFIDENCE SET NOT CALLED BY SVIM


cd ..


echo "Calculating precision and recall is now easy !"


echo "Done"
