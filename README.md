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

- CHROMOSOME: is the chromosome, in the format 'chrN'. Accepted chromosomes are chr1-chr22,chrX,chrY and chrM
- START: where the variant starts
- END: where the variant ends
- ALT: alt type. Possible alt types are 'deletion', 'insertion', 'inversion', 'tr expansion', 'tr contraction', 'ptr', 'atr', 'translocation cut-paste', 'translocation copy-paste' (more details below)
- INFO: info for the alteration (more details below)


## ALT FIELD

VISOR allows users to generate different type of variants specified in the ALT field:

- 'deletion'. Deletes from start(included) to end(included)
- 'insertion'. Inserts a specific sequence immediately after end
- 'inversion'. Inverts from start(included) to end(included)
- 'duplication'. Duplicates from start(included) to end(included)
- 'tr expansion'. Expands an existent tandem repetition. tr is meant to be one of the repetitions in microsatellites regions with START-END pair specified as in the _Examples/GRCh38.microsatellites.bed_ file of this repository (this example is for GRCh38)
- 'tr contraction'. Contracts an existent tandem repetition. Works as described before
- 'ptr'. Inserts a perfect tandem repetition immediately after end
- 'atr' Inserts a approximate tandem repetition immediately after end
- 'translocation cut-paste'. Translocates from start(included) to end(included) to another region. Translocated region is deleted from original position
- 'translocation copy-paste'. Translocates from start(included) to end(included) to another region. Translocated region is not deleted from original position


## INFO FIELD

VISOR allows users to specify parameters in the INFO field:

- INFO for 'deletion' must be _None_
- INFO for 'insertion' must be a valid DNA sequence of any length. Allowed chars are A,C,T,G,N
- INFO for 'inversion' must be _None_
- INFO for 'duplication' must be _number_; number is the number of time segment appears
- INFO for 'tr expansion' must be _motif:number_; motif is a valid DNA motif, number is number of motif to insert
- INFO for 'tr contraction' must be _motif:number_; motif is a valid DNA motif, number is number of motif to delete
- INFO for 'ptr' must be _motif:number_; motif motif is a valid DNA motif, number is number of motif to insert
- INFO for 'atr' must be _motif:number:altnum_; motif is a valid DNA motif, number is number of motif to insert, altnum is the number of alterations; alterations are randomly chosen from 'insertion','deletion','substitution' and each involves one nucleotide only
- INFO for 'translocation cut-paste' must be _haplotype:chromosome:breakpoint:orientation_; haplotype is the haplotype in which region will be translocated ('h1' or 'h2'), chromosome is the chromosome in which region will be translocated (chr1-22, chrX, chrY and chrM are allowed), breakpoint is the number of the base immediately after which translocated region will be put and orientation is the orientation of the sequence ('forward', as it was for the original region, or 'reverse', inverted).
- INFO for 'translocation copy-paste' is the same for 'translocation cut-paste'

### An example .bed file is included in Examples/example.bed

