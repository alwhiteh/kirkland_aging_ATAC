#!/bin/zsh

for i in *.bed; do echo ${i}

findMotifsGenome.pl ${i} dm6 ${i}motifs/ -size 150 -p 8

done 
