# VISOR

![alt text](VISOR.png)


VISOR is a tool written in python for haplotype-specific variants simulations

## Requirements

VISOR requires a unix machine and a working python 3.6 environment. Some python modules are required:

- pyfaidx (v 0.5.5.2)
- pybedtools (v 0.8.0)
- plotly (v 3.7.0)
- pysam (0.15.2) - for the sscounter.py script -

These modules will be automatically installed during VISOR setup and nothing from the user is required.

Moreover, for simulations, VISOR requires working installations of:

- samtools (https://github.com/samtools/samtools)
- wgsim (https://github.com/lh3/wgsim)
- pbsim (https://github.com/pfaucon/PBSIM-PacBio-Simulator)
- bwa (https://github.com/lh3/bwa)
- minimap2 (https://github.com/lh3/minimap2)

If you do not have a py36 environment, it can be generated with anaconda (https://docs.anaconda.com/anaconda):

```sh 

conda create -n py36 python=3.6 anaconda

```

If the aforementioned tools are not installed in your path, they can all be installed with anaconda:


```sh 

source activate py36
conda install -c bioconda samtools 
conda install -c bioconda wgsim
conda install -c bioconda pbsim
conda install -c bioconda bwa
conda install -c bioconda minimap2 
#or in one go:
# conda install -c bioconda samtools wgsim pbsim bwa minimap2
```

## Install VISOR

```sh
source activate py36 #or conda activate py36 for more recent versions
git clone https://github.com/davidebolo1993/VISOR.git
cd VISOR
python setup.py install

```

## VISOR Submodules

VISOR is built on 3 submodules:

- __VISOR HACk__: generates .fa haplotypes with user-defined variants
- __VISOR SHORtS__: simulate double-strand or single-strand short-reads .bam files with variants
- __VISOR LASeR__: simulate long-reads .bam files with variants


## VISOR HACk

```sh
VISOR HACk -h #print help

VISOR HACk -g genome.fa -bed bed1.bed bed2.bed -o hackout

```

Inputs to VISOR HACk are:

- genome.fa is the reference genome in .fasta format
- bed1.bed is and bed2.bed (can be also more) are .bed files containing variants to create for each haplotype.


VISOR HACk outputs a fasta (.fa) file with specified variants for each .bed in the output folder (hackout/h1/h1.fa and hackout/h2/h2.fa in this example) 


#### Construct the .bed file/s

.bed file/s must contain 5 columns WITHOUT header: __CHROMOSOME__, __START__, __END__, __ALT__, __INFO__

- __CHROMOSOME__: is the chromosome, in the format 'chrN'. Accepted chromosomes are the ones also present in the reference genome
- __START__: where the variant starts
- __END__: where the variant ends
- __ALT__: alt type. Possible alt types are 'deletion', 'insertion', 'inversion', 'tandem duplication', 'SNP', 'tandem repeat expansion', 'tandem repeat contraction', 'perfect tandem repetition', 'approximate tandem repetition', 'translocation cut-paste', 'translocation copy-paste' (or 'interspersed duplication'), 'reciprocal translocation' (more details below)
- __INFO__: informations for the alteration (more details below)


##### ALT FIELD

VISOR HACk allows users to generate different type of variants specified in the ALT field of the .bed files:

- __'deletion'__. Deletes from start (included) to end (included)
- __'insertion'__. Inserts a specific sequence immediately after end
- __'inversion'__. Inverts from start (included) to end (included)
- __'tandem duplication'__. Duplicates from start (included) to end (included) immediately after end
- __'SNP'__. Introduces a single nucleotide polimorphism in end
- __'tandem repeat expansion'__. Expands an existent tandem repetition. tandem repeat is meant to be one of the repetitions present in microsatellites regions with START-END pair specified as in the _Examples/GRCh38.microsatellites.bed_ file of this repository (this example is taken from UCSC for GRCh38 reference)
- __'tandem repeat contraction'__. Contracts an existent tandem repetition. Works as described for 'tandem repeat expansion'
- __'perfect tandem repetition'__. Inserts a perfect tandem repetition immediately after end
- __'approximate tandem repetition'__ Inserts a approximate tandem repetition immediately after end
- __'translocation cut-paste'__. Translocates from start (included) to end (included) to another region. Translocated region is deleted from original position
- __'translocation copy-paste'__ or __'interspersed duplication'__. Translocates from start (included) to end (included) to another region. Translocated region is not deleted from original position
- __'reciprocal translocation'__. Translocates from start (included) to end (included) to another region and translocates the destination region back to the first one.


##### INFO FIELD

VISOR HACk requires some users-defined parameteres in the INFO field of the .bed files:

- INFO for __'deletion'__ must be __None__
- INFO for __'insertion'__ must be a valid DNA sequence of any length. Allowed chars are A,C,T,G,N
- INFO for __'inversion'__ must be __None__
- INFO for __'tandem duplication'__ must be __number__; number is the number of time segment appears
- INFO for __'SNP'__ must be __nuc__; nuc is the nucleotide that will be used to introduce the variant
- INFO for __'tandem repeat expansion'__ must be __motif:number__ motif is a valid DNA motif, number is number of motif to insert
- INFO for __'tandem repeat contraction'__ must be __motif:number__; motif is a valid DNA motif, number is number of motif to delete
- INFO for __'perfect tandem repetition'__ must be __motif:number__; motif motif is a valid DNA motif, number is number of motif to insert
- INFO for __'approximate tandem repetition'__ must be __motif:number:alterations__ ; motif is a valid DNA motif, number is number of motif to insert, alterations is the number of alterations; alterations are randomly chosen from 'insertion','deletion','substitution' and each involves one nucleotide only
- INFO for __'translocation cut-paste'__ must be __haplotype:chromosome:breakpoint:orientation__; haplotype is the haplotype in which region will be translocated ('h1', 'h2', ...), chromosome is the chromosome in which region will be translocated (chr1-22, chrX, chrY and chrM are allowed), breakpoint is the number of the base immediately before the one where translocated region will start and orientation is the orientation of the sequence ('forward', if the orientation should be the same of the original region, or 'reverse', if the orientation should be inverted).
- INFO for __'translocation copy-paste'__ is the __same for 'translocation cut-paste'__
- INFO for __'reciprocal translocation'__ is the __haplotype:chromosome:breakpoint:orientation1:orientation2__; haplotype is the haplotype in which region will be translocated ('h1', 'h2', ...), chromosome is the chromosome in which region will be translocated; breakpoint is the number of the base immediately before the one where translocated region will start; orientation1 is the orientation of the first region ('forward', if the orientation should be the same of the original region, or 'reverse', if the orientation should be inverted) and orientation2 is the orientation of the second region.


## VISOR SHORtS and VISOR LASeR

```sh

VISOR SHORtS -h #print help

VISOR SHORtS -g genome.fa -hapfa h1.fa h2.fa -bed sim.bed -o shortsout #default double-strand sequencing simulations

VISOR LASeR -h #print help

VISOR LASeR -g genome.fa -hapfa h1.fa h2.fa -bed sim.bed -o laserout #long-reads sequencing simulations

```

Inputs to VISOR SHORtS and VISOR LASeR are:

- genome.fa is the reference genome in .fasta format
- h1.fa and h2.fa (can be also more) are the haplotype .fasta file generated with VISOR HACk
- sim.bed is the .bed file containing regions to simulate from h1.fa and h2.fa like the one in _Examples/SHORtS.LASeR.bed_. If multiple regions are given, the .bam files generated for each region will then be merged.

Additionally, both VISOR SHORtS and VISOR LASeR allow to specify an allelic fraction (-af), which is the fraction of reads that for the simulated region in sim.bed files will carry the structural variants.

.bed file must contain 4 columns without header: __CHROMOSOME__, __START__, __END__, __COVERAGE BIAS__

- __CHROMOSOME__: is the chromosome, in the format 'chrN'. Accepted chromosomes are chr1-chr22, chrX, chrY and chrM
- __START__: start position for the region that will be simulated
- __END__: end position for the region that will be simulated
- __COVERAGE BIAS__: a float to specify a deviation from the wanted coverage (80.0 means that the region is covered by the 80% of the reads that were supposed to cover the region).

VISOR SHORtS and VISOR LASeR outputs a .srt.bam file (and its index) for each region in the haplotype-specific .bed files. 

#### VISOR SHORtS and VISOR LASeR Simulations

Simulations for short-reads data are run using wgsim, most of which parameters can be specified by the user, bwa and samtools. Simulations for long-reads data are run using pbsim, most of which parameters can be specified by the user (the model_qc_clr file required is included in VISOR), minimap2 and samtools.


## VISOR SHORtS for single-strand (strand-seq) simulations

VISOR SHORtS can simulate single-strand (strand-seq) .bam files

```sh

VISOR SHORtS -h #print help

VISOR SHORtS -g genome.fa -hapfa h1.fa h2.fa -bed sim.bed -t single-strand -c 1 -o shortsout #single-strand simulations without noise

VISOR SHORtS -g genome.fa -hapfa h1.fa h2.fa -bed sim.bed -t single-strand -c 1 -n 5 -o shortsout #single-strand simulations with 5 % of read pairs with incorrect orientation

VISOR SHORtS -g genome.fa -hapfa h1.fa h2.fa -bed sim.bed -t single-strand -c 1 -n 5 -scebed sce.bed -o shortsout #single-strand simulations with 5 % of read pairs with incorrect orientation and sister chromatid exchange for regions and haplotypes specified in sce.bed


```

When working in single-strand mode (-t single-strand), for each region specified in the .bed file, VISOR outputs 2 .bam files: in the 'watson' .bam file, read 1 and read 2 pairs have all forward and reverse orientation respectively; in the 'crick' .bam file, read 1 and read 2 pairs have all reverse and forward orientation respectively. It is also possible to specify a percentage of noise (pairs with incorrect orientation) that will be included in the crick and watson .bam files using the -n parameter. Users can also specify a .bed files (-scebed) containing informations (chromosome, start, end, haplotype) for haplotypes and regions in which sister chromatid exchange will be performed.

#### Plot read pairs count

When simulating haplotype-specific variants in single-strand simulations, one can re-create the possible inherited template strands in daughter cells by combining the watson or crick .bam file for 1 haplotype (W1 or C1) with  the watson or crick  one from the other(W2 or C2), for example by using samtools merge

- W1 and W2
- W1 and C2
- C1 and W2
- C1 and C2

VISOR/scripts/sscounter.py offers the possibility to plot an interactive visualization of the read pairs count for given .bam file/s

```sh

python VISOR/scripts/sscounter.py -h #print help

python VISOR/scripts/sscounter.py -g genome.fa -bam .bam1 .bam2  -O sscounterout #generate an .html that, for each chromosome,  compares read pairs count for watson and crick strands of the given .bam files

```
