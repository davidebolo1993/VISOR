# VISOR
VarIant SimulatOR


VISOR is a small tool written in python for fast haplotype-specific variants simulation.

## Requirements

VISOR requires a working python 3.6 environment and depends on the following python packages:

- pyfaidx (v 0.5.5.2)
- pybedtools (v 0.8.0)

## Install VISOR

```sh
git clone https://github.com/davidebolo1993/VISOR.git
cd VISOR
python setup.py install

```

## Run VISOR

```sh

VISOR -g genome.fa -bedh1 bedh1.bed -bedh2 bedh2.bed -O pathout

```


## Inputs

Inputs to VISOR are:

- a genome .fasta file
- a .bed file containing variants for haplotype 1
- a .bed file containing variants for haplotype 2 (optional)

.bed file must contain 5 columns without header: CHROMOSOME, START, END, ALT, INFO

- __CHROMOSOME__: is the chromosome, in the format 'chrN'. Accepted chromosomes are chr1-chr22,chrX,chrY and chrM
- __START__: where the variant starts
- __END__: where the variant ends
- __ALT__: alt type. Possible alt types are 'deletion', 'insertion', 'inversion', 'duplication', 'snp', 'tr expansion', 'tr contraction', 'ptr', 'atr', 'translocation cut-paste', 'translocation copy-paste' (more details below)
- __INFO__: info for the alteration (more details below)

_An example .bed file is included in Examples/example.bed_


## Outputs

A fasta (.fa) for each haplotype in the output folder (pathout/h1.fa and pathout/h2.fa) containing specified alterations.

## ALT FIELD

VISOR allows users to generate different type of variants specified in the ALT field:

- __'deletion'__. Deletes from start (included) to end (included)
- __'insertion'__. Inserts a specific sequence immediately after end
- __'inversion'__. Inverts from start (included) to end (included)
- __'duplication'__. Duplicates from start (included) to end (included) immediately after end
- __'snp'__. Introduces a single nucleotide polimorphism in end
- __'tr expansion'__. Expands an existent tandem repetition. tr is meant to be one of the repetitions in microsatellites regions with START-END pair specified as in the _Examples/GRCh38.microsatellites.bed_ file of this repository (this example is for GRCh38)
- __'tr contraction'__. Contracts an existent tandem repetition. Works as described before
- __'ptr'__. Inserts a perfect tandem repetition immediately after end
- __'atr'__ Inserts a approximate tandem repetition immediately after end
- __'translocation cut-paste'__. Translocates from start (included) to end (included) to another region. Translocated region is deleted from original position
- __'translocation copy-paste'__. Translocates from start (included) to end (included) to another region. Translocated region is not deleted from original position


## INFO FIELD

VISOR requires some users-defined parameteres in the INFO field:

- INFO for 'deletion' must be __None__
- INFO for 'insertion' must be a valid DNA sequence of any length. Allowed chars are A,C,T,G,N
- INFO for 'inversion' must be __None__
- INFO for 'duplication' must be __number__; number is the number of time segment appears
- INFO for 'snp' must be __nuc__; nuc is the nucleotide that will be used to introduce the variant
- INFO for 'tr expansion' must be __motif:number__ motif is a valid DNA motif, number is number of motif to insert
- INFO for 'tr contraction' must be __motif:number__; motif is a valid DNA motif, number is number of motif to delete
- INFO for 'ptr' must be __motif:number__; motif motif is a valid DNA motif, number is number of motif to insert
- INFO for 'atr' must be __motif:number:altnum__; motif is a valid DNA motif, number is number of motif to insert, altnum is the number of alterations; alterations are randomly chosen from 'insertion','deletion','substitution' and each involves one nucleotide only
- INFO for 'translocation cut-paste' must be __haplotype:chromosome:breakpoint:orientation__; haplotype is the haplotype in which region will be translocated ('h1' or 'h2'), chromosome is the chromosome in which region will be translocated (chr1-22, chrX, chrY and chrM are allowed), breakpoint is the number of the base immediately after which translocated region will be put and orientation is the orientation of the sequence ('forward', if the orientation should be the same of the original region, or 'reverse', if the orientation should be inverted).
- INFO for 'translocation copy-paste' is the same for 'translocation cut-paste'
