#!/bin/bash

cd $1

for f in *fastq; do

  sed -i -e's/\/[0-9]//g' $f && gzip $f

done
