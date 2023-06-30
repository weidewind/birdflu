#!/usr/bin/bash

full='false'
while getopts i:f: flag
do
    case "${flag}" in
        i) cdhit_id=${OPTARG};;
		f) full='true';;
    esac
done

ci=$(echo $cdhit_id | sed -r 's/\.//g')
datafolder=/export/home/popova/workspace/birdflu/data
outputfolder=/export/home/popova/workspace/birdflu/output
alifolder=${outputfolder}/aligned_mafft_cdhit_${ci}_indexing
scriptfolder=/export/home/popova/workspace/birdflu/scripts
metapath=${datafolder}/meta_short.csv
blastfolder=/export/home/popova/workspace/birdflu/output/blast
remove=/export/home/popova/workspace/birdflu/data/remove.list
nsubs=(N4 N5 N1 N2 N3 N6 N7 N8 N9)

cd $datafolder
mkdir $alifolder


nuclfasta=nucl.slow10it.edited.${ci}.cdhit.realigned
# protfasta=prot.slow10it.edited.${ci}.cdhit.realigned
# nuclconcat=${alifolder}/concat.${nuclfasta}

prep_and_pal2nal() {
	nucl=$1
	msa=$2
	output=$3
	## pal2nal does not understand Ns (but does understand dots)
	perl -pi -e 's/N/./g unless $_ =~ "^>.*"' $1
	## pal2nal does not understand gapped codons (eg -TG), thus leading and trailing gaps in, respectively, leading and trailing codons are changed to dots
	## since we aligned seqs with reference proteins by codon-aware virulign, we do not expected any internal codons with gaps
	python ${scriptfolder}/fix_gapfices.py --input $1 --output ${1}.dotted
	pal2nal.pl $msa ${1}.dotted -output fasta >${output} 2>${output}.err
	unlink ${1}.dotted
}

mafft_and_pal2nal() {
	## nucl fasta
	nfasta=$1
	## output file name
	nmsa=$2
	it=$3
	python ${scriptfolder}/transeq.py --input $nfasta --output $nfasta.transeq
	mafft --anysymbol --thread -1 --maxiterate $it --globalpair $nfasta.transeq >$nfasta.transeq.slow${it}it
	prep_and_pal2nal $nfasta $nfasta.transeq.slow${it}it $nmsa
	}

## Merge intra-subtype alignments (without seed) 
mafft_merge_and_pal2nal() {
	local nuclfasta=$1 # extension
	shift
	nsubs=("$@")
	
	protconcat=${alifolder}/concat.${nuclfasta}.transeq
	nuclconcat=${alifolder}/concat.${nuclfasta}
	>$protconcat
	>$nuclconcat
	>${alifolder}/prot.subMSAtable
	seqcount=1
	for nsub in "${nsubs[@]}"
		do 
			python ${scriptfolder}/transeq.py --input ${alifolder}/${nsub}.${nuclfasta} --output ${alifolder}/${nsub}.${nuclfasta}.transeq
			cat ${alifolder}/${nsub}.${nuclfasta}.transeq >>$protconcat
			cat ${alifolder}/${nsub}.${nuclfasta} >>$nuclconcat
			n=$(grep -c ">" ${alifolder}/${nsub}.${nuclfasta}.transeq)
			t=$(($n+$seqcount))
			echo $(seq $seqcount $(($t-1)) | tr '\n' ' ' | sed 's/\s$//g') >>${alifolder}/prot.subMSAtable
			seqcount=$t
		done
	## seems that --maxiterate does not work without --globalpair (or --localpair)
	## without --anysymbol mafft changes * to a gap
	## refinement is indispensable, since our sampling is very uneven (>50% seqs are from zaire ebolavirus) 
	echo "Merging alignments for subtypes with mafft+pal2nal.."
	mafft --anysymbol --thread -1 --maxiterate 50 --globalpair --merge ${alifolder}/prot.subMSAtable $protconcat >${protconcat}.msa
	## We already prepared dotted file for pal2nal
	pal2nal.pl ${protconcat}.msa $nuclconcat -output fasta >${nuclconcat}.msa 2>${nuclconcat}.msa.err

}

## Merge intra-subtype alignments (with seed) 
mafft_merge_with_seed_and_pal2nal() {
	local nuclfasta=$1 # extension
	local seed=$2
	local it=$3
	shift
	nsubs=("$@")
	
	protconcat=${alifolder}/concat.${nuclfasta}.transeq
	nuclconcat=${alifolder}/concat.${nuclfasta}
	>$protconcat
	>$nuclconcat
	>${alifolder}/prot.withseed.subMSAtable
	seqcount=$((1+$(grep -c ">" $seed)))
	for nsub in "${nsubs[@]}"
		do 
			python ${scriptfolder}/transeq.py --input ${alifolder}/${nsub}.${nuclfasta} --output ${alifolder}/${nsub}.${nuclfasta}.transeq
			cat ${alifolder}/${nsub}.${nuclfasta}.transeq >>$protconcat
			cat ${alifolder}/${nsub}.${nuclfasta} >>$nuclconcat
			n=$(grep -c ">" ${alifolder}/${nsub}.${nuclfasta}.transeq)
			t=$(($n+$seqcount))
			echo $(seq $seqcount $(($t-1)) | tr '\n' ' ' | sed 's/\s$//g') >>${alifolder}/prot.withseed.subMSAtable
			seqcount=$t
		done
	## seems that --maxiterate does not work without --globalpair (or --localpair)
	## without --anysymbol mafft changes * to a gap
	## refinement is indispensable, since our sampling is very uneven (>50% seqs are from zaire ebolavirus) 
	echo "Merging alignments for subtypes with mafft+pal2nal.."
	mafft --anysymbol --thread -1 --seed $seed --maxiterate $it --globalpair --merge ${alifolder}/prot.subMSAtable $protconcat >${protconcat}.withseed.msa
	## We already prepared dotted file for pal2nal
	pal2nal.pl ${protconcat}.withseed.msa $nuclconcat -output fasta >${nuclconcat}.withseed.msa 2>${nuclconcat}.withseed.msa.err
}


# ## take aligned seqs (now you can just take fasta, because we filter it later on), throw out gaps, filter the resulting fasta, cd-hit and mafft 
# for nsub in "${nsubs[@]}"
	# do
	# ## Identity 85%, local alignment, 100% coverage for the shorter seq, any coverage for the longer seq, 20 threads
	# ## cdhit cannot work with gaps(
	# ## TODO make seq filtering at the early stages - and then you won't need to sed and redo the analysis
	# ## ${alifolder}/${nsub}.nucl.slow10it.edited is an output from cdhit.sh, after manual correction. This could be just a clean fasta!
	# ## TODO filtering and clustering and deleteing seqs which cluster with some other subtype\species (misspecified seqs)
	# cp /export/home/popova/workspace/birdflu/output/aligned_mafft_cdhit/${nsub}.nucl.slow10it.edited ${alifolder}/${nsub}.nucl.slow10it.edited
	# perl -p -e 's/-//g unless $_ =~ "^>.*"' ${alifolder}/${nsub}.nucl.slow10it.edited >${alifolder}/${nsub}.nucl.slow10it.edited.nogaps.unfiltered
	# sed -i "/^$/d" ${alifolder}/${nsub}.nucl.slow10it.edited.nogaps.unfiltered

	# faFilter -namePatList=$remove -v ${alifolder}/${nsub}.nucl.slow10it.edited.nogaps.unfiltered ${alifolder}/${nsub}.nucl.slow10it.edited.nogaps
	# cd-hit -c $cdhit_id -G 0 -aS 1 -aL 0 -T 20 -i ${alifolder}/${nsub}.nucl.slow10it.edited.nogaps -o ${alifolder}/${nsub}.nucl.slow10it.edited.${ci}.cdhit
	# mafft_and_pal2nal ${alifolder}/${nsub}.nucl.slow10it.edited.${ci}.cdhit ${alifolder}/${nsub}.${nuclfasta} 50
	# done


# ## Merge cdhitted intra-subtype alignments (without seed) 
# ## nuclfasta - ext for cdhitted intra-subtype alignments
# ## Output: ${alifolder}/concat.${nuclfasta}.msa
# mafft_merge_and_pal2nal $nuclfasta "${nsubs[@]}"

# ## Merge full intra-subtype alignments (with cdhitted alignment as a seed) 
# ## nuclfasta - ext for full intra-subtype alignments
# if ${full}
	# then
	# mafft_merge_with_seed_and_pal2nal nucl.slow10it.edited ${alifolder}/concat.${nuclfasta}.msa 10 "${nsubs[@]}"
	# fi


# ## Make consensus for eahc subtype
>${outputfolder}/subtype_alignments/all.consensus.fasta
for nsub in "${nsubs[@]}"
	do
	# cons ${outputfolder}/subtype_alignments/${nsub}.nucl.slow10it.filtered ${outputfolder}/subtype_alignments/${nsub}.consensus
	# echo "Py consensus for $nsub.."
	# python ${scriptfolder}/consensus.py --msa ${outputfolder}/subtype_alignments/${nsub}.nucl.slow10it.filtered --name ${nsub}_consensus --output ${outputfolder}/subtype_alignments/${nsub}.py.consensus
	# sed -i 's/-//g' ${outputfolder}/subtype_alignments/${nsub}.py.consensus
	cat ${outputfolder}/subtype_alignments/${nsub}.py.consensus >>${outputfolder}/subtype_alignments/all.consensus.fasta
	done


numid=OQ632858_NA # Kamil's seq
window=15
out_cons_thr=0.2
pivot_subtype=N1
## Search for windows where N1 is conservative and all other groups have small frequency of N1 consensus symbol
## divergence file: 
## 1. index according to numid (last char of the window)
## 2. mean consensus char freq in N1 across window positions
## 3. mean N1 consensus char freq across all other groups and all window positions
## .wg modification means that gaps are treated the same as other symbols
## Add inexing sequence 
# mafft --add refs/indexing_seq ${alifolder}/concat.${nuclfasta}.msa >${alifolder}/concat.${nuclfasta}.msa.temp

## Add consensus sequence
for nsub in N2 #"${nsubs[@]}"
	do
	mafft --add ${outputfolder}/subtype_alignments/all.consensus.fasta ${alifolder}/concat.${nuclfasta}.msa >${alifolder}/concat.${nuclfasta}.msa.withcons
	done

# ## Need to filter for cdhit.1, because it was produced by another script 
# faFilter -namePatList=$remove -v ${alifolder}/concat.${nuclfasta}.msa.temp ${alifolder}/concat.${nuclfasta}.msa
# python ${scriptfolder}/write_divergence_table.py --pivot $pivot_subtype --msa ${alifolder}/concat.${nuclfasta}.msa --meta ${datafolder}/meta_short.csv \
# --output ${alifolder}/${pivot_subtype}.divergence.table.wg.${numid} --numeration_id $numid

# ## Divergence from all other subtypes together
# python ${scriptfolder}/search_divergent.py --table ${alifolder}/${pivot_subtype}.divergence.table.wg.${numid} \
# --window $window --out_cons_thr $out_cons_thr --output ${alifolder}/${pivot_subtype}.divergence.wg.window.${window}.outthr${out_cons_thr}.${numid}
# ## Divergence from all other subtypes individually
# for nsub in "${nsubs[@]}"
	# do
	# python ${scriptfolder}/search_divergent.py --table ${alifolder}/${pivot_subtype}.divergence.table.wg.${numid} \
	# --out_cons_thr $out_cons_thr  --window $window --groups $nsub \
	# --output ${alifolder}/${pivot_subtype}.divergence.wg.window.${window}.outthr${out_cons_thr}.${numid}.comparedto.${nsub}
	# done


cd $datafolder

## Check how the MANUALLY selected region maps to other sequences

## Create blastdb
# ./${scriptfolder}/filter.sh -f $datafolder -i 1.cdhit.fasta -o 1.cdhit.filtered.fasta
# for nsub in "${nsubs[@]}"
	# do
	# mv ${datafolder}/${nsub}.1.cdhit.filtered.fasta ${blastdb}/${nsub}.1.cdhit.filtered.fasta
	# done
# cd $BLASTDB
# for nsub in "${nsubs[@]}"
	# do
	# makeblastdb -in ${datafolder}/subtype_fasta/${nsub}.1.cdhit.filtered.fasta -parse_seqids -blastdb_version 5 \
	# -title ${nsub} -out ${nsub} -dbtype nucl
# done


#loc=$(echo ${start}-${end}) #191-206
mkdir ${alifolder}/mafft_check
mkdir ${alifolder}/seqkit_check
mkdir ${alifolder}/bowtie_check
mkdir ${alifolder}/bowtie_check_minscore_letter

bowtie2_index=/export/home/popova/local/bowtie2_index

## Create bowtie2 index
# cd $bowtie2_index
# for nsub in "${nsubs[@]}"
	# do
	# bowtie2-build -f ${datafolder}/subtype_fasta/${nsub}.1.cdhit.filtered.fasta $nsub
	# done
# cd $datafolder

# locs=(1000-1022 1141-1155 1177-1201 1212-1253 130-145 1363-1388 226-247 259-275 339-347 378-398 726-750)

# for loc in "${locs[@]}"
	# do
		# start=$(echo $loc | cut -d"-" -f1)
		# end=$(echo $loc | cut -d"-" -f2)
		# mm=$(((${end}-${start})/5))
		# echo $loc
		# echo $mm
		# >${alifolder}/bowtie_check_minscore_letter/${numid}.${loc}.matchcount
		
		# ## Cut selected region
		# python ${scriptfolder}/cut_fragment.py --fasta refs/${numid}.fasta --start $start --stop $end --output refs/${numid}.${loc}.fasta

		# ## Map selected regions to other seqs
		# for nsub in "${nsubs[@]}"
			# do
			# ## --score-min: defines score function parameters, for L,-20.0,-1.0 min_score= -20-1*L (L=algn length), which allows ~6 mismatches (mm penalty = 6)
			# ## -L seed length; -N 1: 1 mismatch allowed in the seed
			# bowtie2 --score-min L,0,-1.2 -L 7 -N 1 -a -x ${bowtie2_index}/$nsub -f refs/${numid}.${loc}.fasta >${alifolder}/bowtie_check_minscore_letter/${nsub}.${numid}.${loc}.out #--suppress 1,5,6,7
			# echo $nsub >>${alifolder}/bowtie_check_minscore_letter/${numid}.${loc}.matchcount
			# grep -o -E "NM:i:[0-9]+" ${alifolder}/bowtie_check_minscore_letter/${nsub}.${numid}.${loc}.out | sort | uniq -c >>${alifolder}/bowtie_check_minscore_letter/${numid}.${loc}.matchcount
			# # seqkit locate -m $mm -i -f refs/${numid}.${loc}.fasta ${datafolder}/subtype_fasta/${nsub}.1.cdhit.filtered.fasta >${alifolder}/seqkit_check/${nsub}.${numid}.${loc}.out
			# # mafft --addfragments refs/${numid}.${loc}.fasta ${outputfolder}/subtype_alignments/${nsub}.nucl.slow10it.filtered >${alifolder}/mafft_check/${nsub}.${numid}.${loc}.out
			# # blastn -task blastn-short -db $nsub -query refs/${numid}.fasta -query_loc $loc \
			# # -out ${alifolder}/blast_check/${nsub}.${numid}.${loc}.out
			# # bwa index -a is ${datafolder}/subtype_fasta/${nsub}.1.cdhit.filtered.fasta -p ${nsub}
			# # bwa aln -n $mm ${nsub} refs/${numid}.${loc}.fasta >${alifolder}/bwa_check/${nsub}.${numid}.${loc}.sai
			# # bwa samse ${nsub} ${alifolder}/bwa_check/${nsub}.${numid}.${loc}.sai refs/${numid}.${loc}.fasta >${alifolder}/bwa_check/${nsub}.${numid}.${loc}.sam
			# done
	# done
