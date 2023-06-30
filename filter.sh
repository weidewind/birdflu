#!/usr/bin/bash

while getopts i:o:f: flag
do
    case "${flag}" in
        i) in_suff=${OPTARG};;
		o) out_suff=${OPTARG};;
		f) folder=${OPTARG};;
    esac
done

remove=/export/home/popova/workspace/birdflu/data/remove.list
nsubs=(N4 N5 N1 N2 N3 N6 N7 N8 N9)
cd $folder

## take aligned seqs (now you can just take fasta, because we filter it later on), throw out gaps, filter the resulting fasta, cd-hit and mafft 
for nsub in "${nsubs[@]}"
	do
	faFilter -namePatList=$remove -v ${nsub}.${in_suff} ${nsub}.${out_suff}
	done

