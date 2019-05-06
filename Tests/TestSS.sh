echo "Creating test files"
mkdir SStest && cd SStest

echo "Download reference"
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa chr22 > chr22.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa chr11 > chr11.fa
cat chr22.fa chr11.fa > chr11_22.fa


echo "Download scripts from VISOR git"


wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/scripts/ssmerger.py
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/scripts/sscounter.py



echo "Getting files for simulations from VISOR git"


mkdir files && cd files

mkdir hetdup && cd hetdup

wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/hetdup/VISOR.h1.SVs.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/hetdup/VISOR.h2.SVs.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/hetdup/VISOR.sim.bed


cd ..

mkdir hetdel && cd hetdel

wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/hetdel/VISOR.h1.SVs.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/hetdel/VISOR.h2.SVs.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/hetdel/VISOR.sim.bed


cd ..


mkdir hetinv && cd hetinv

wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/hetinv/VISOR.h1.SVs.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/hetinv/VISOR.h2.SVs.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/hetinv/VISOR.sim.bed


cd ..



mkdir hetinvdup && cd hetinvdup


wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/hetinvdup/VISOR.h1.SVs.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/hetinvdup/VISOR.h2.SVs.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/hetinvdup/VISOR.sim.bed


cd ..


mkdir transloc && cd transloc

wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/transloc/VISOR.h1.SVs.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/transloc/VISOR.h2.SVs.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/transloc/VISOR.sim.bed



cd ..


mkdir SCE && cd SCE

wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/SCE/VISOR.h1.SVs.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/SCE/VISOR.h2.SVs.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/SCE/VISOR.sim.bed
wget https://raw.githubusercontent.com/davidebolo1993/VISOR/master/Tests/TestSS/SCE/VISOR.h1.SCE.bed


cd ..


cd ..


echo "Done"


echo "Simulating strand-seq data with 5% noise for TANDEM DUPLICATION"

mkdir hetdup

VISOR HACk -g chr22.fa -bed files/hetdup/VISOR.h1.SVs.bed files/hetdup/VISOR.h2.SVs.bed -o hetdup/Templates
VISOR SHORtS -g chr22.fa -s hetdup/Templates -bed files/hetdup/VISOR.sim.bed -c 1 -t single-strand -n 5.0 -o hetdup/Results
python ssmerger.py -f hetdup/Results/0 hetdup/Results/1 -s W C -o hetdup/Final
python sscounter.py -g GRCh38_full_analysis_set_plus_decoy_hla.fa -bam hetdup/Final/WC.srt.bam -o hetdup/Final -l H1_W.H2_C.H2_DUP -c chr22


echo "Done"

echo "Simulating strand-seq data with 5% noise for DELETION"

mkdir hetdel

VISOR HACk -g chr22.fa -bed files/hetdel/VISOR.h1.SVs.bed files/hetdel/VISOR.h2.SVs.bed -o hetdel/Templates
VISOR SHORtS -g chr22.fa -s hetdel/Templates -bed files/hetdel/VISOR.sim.bed -c 1 -t single-strand -n 5.0 -o hetdel/Results
python ssmerger.py -f hetdel/Results/0 hetdel/Results/1 -s W C -o hetdel/Final
python sscounter.py -g GRCh38_full_analysis_set_plus_decoy_hla.fa -bam hetdel/Final/WC.srt.bam -o hetdel/Final -l H1_W.H2_C.H2_DEL -c chr22

echo "Done"


echo "Simulating strand-seq data with 5% noise for INVERSION"

mkdir hetinv

VISOR HACk -g chr22.fa -bed files/hetinv/VISOR.h1.SVs.bed files/hetinv/VISOR.h2.SVs.bed -o hetinv/Templates
VISOR SHORtS -g chr22.fa -s hetinv/Templates -bed files/hetinv/VISOR.sim.bed -c 1 -t single-strand -n 5.0 -o hetinv/Results
python ssmerger.py -f hetinv/Results/0 hetinv/Results/1 -s W C -o hetinv/Final
python sscounter.py -g GRCh38_full_analysis_set_plus_decoy_hla.fa -bam hetinv/Final/WC.srt.bam -o hetinv/Final -l H1_W.H2_C.H2_INV -c chr22


echo "Done"

echo "Simulating strand-seq data with 5% noise for INVERTED DUPLICATION"

mkdir hetinvdup

VISOR HACk -g chr22.fa -bed files/hetinvdup/VISOR.h1.SVs.bed files/hetinvdup/VISOR.h2.SVs.bed -o hetinvdup/Templates
VISOR SHORtS -g chr22.fa -s hetinvdup/Templates -bed files/hetinvdup/VISOR.sim.bed -c 1 -t single-strand -n 5.0 -o hetinvdup/Results
python ssmerger.py -f hetinvdup/Results/0 hetinvdup/Results/1 -s W C -o hetinvdup/Final
python sscounter.py -g GRCh38_full_analysis_set_plus_decoy_hla.fa -bam hetinvdup/Final/WC.srt.bam -o hetinvdup/Final -l H1_W.H2_C.H2_INVDUP -c chr22


echo "Done"

echo "Simulating strand-seq data with 5% noise for BALANCED TRANSLOCATION"

mkdir transloc

VISOR HACk -g chr11_22.fa -bed files/transloc/VISOR.h1.SVs.bed files/transloc/VISOR.h2.SVs.bed -o transloc/Templates
VISOR SHORtS -g chr11_22.fa -s transloc/Templates -bed files/transloc/VISOR.sim.bed -c 1 -t single-strand -n 5.0 -o transloc/Results
python ssmerger.py -f transloc/Results/0 transloc/Results/1 -s W C -o transloc/Final
python sscounter.py -g GRCh38_full_analysis_set_plus_decoy_hla.fa -bam transloc/Final/WC.srt.bam -o transloc/Final -l H1_W.H2_C.22_11_TRANSLOC -c chr22 chr11


echo "Done"


echo "Simulating strand-seq data with 5% and SCE in haplotype 1"

mkdir SCE

VISOR HACk -g chr22.fa -bed files/SCE/VISOR.h1.SVs.bed files/SCE/VISOR.h2.SVs.bed -o SCE/Templates
VISOR SHORtS -g chr22.fa -s SCE/Templates -bed files/SCE/VISOR.sim.bed -c 1 -t single-strand -n 5.0 -scebed files/SCE/VISOR.h1.SCE.bed -o SCE/Results
python ssmerger.py -f SCE/Results/0 SCE/Results/1 -s W C -o SCE/Final
python ssmerger.py -f SCE/Results/0 SCE/Results/1 -s C W -o SCE/Final
python sscounter.py -g GRCh38_full_analysis_set_plus_decoy_hla.fa -bam SCE/Final/WC.srt.bam -o SCE/Final/A -l H1_W.H2_C.H1_SCE -c chr22
python sscounter.py -g GRCh38_full_analysis_set_plus_decoy_hla.fa -bam SCE/Final/CW.srt.bam -o SCE/Final/B -l H1_C.H2_W.H1_SCE -c chr22

echo "Done"

