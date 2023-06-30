#!/usr/bin/bash
datafolder=/export/home/popova/workspace/birdflu/data
outfolder=/export/home/popova/workspace/birdflu/output
scriptfolder=/export/home/popova/workspace/birdflu/scripts
metapath=${datafolder}/meta_short.csv
cd $datafolder
nsubs=(N4 N5 N1 N2 N3 N6 N7 N8 N9)



## Select fasta for each N subtype
## Virulign needs path export: screen -U -S virnp env LD_LIBRARY_PATH=$LD_LIBRARY_PATH 
for nsub in "${nsubs[@]}"
	do
	mafft --globalpair --maxiterate 10 --thread -1 --adjustdirection ${nsub}.1.cdhit.fasta >${outfolder}/aligned/${nsub}.slow.10it.msa 2>${outfolder}/aligned/${nsub}.slow.10it.msa.err
	done
