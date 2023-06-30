#!/usr/bin/bash

datafolder=/export/home/popova/workspace/birdflu/data
alifolder=/export/home/popova/workspace/birdflu/output/aligned_2
scriptfolder=/export/home/popova/workspace/birdflu/scripts
metapath=${datafolder}/meta_short.csv
blastfolder=/export/home/popova/workspace/birdflu/output/blast
cd $datafolder
nsubs=(N4 N5 N1 N2 N3 N6 N7 N8 N9)


# for nsub in "${nsubs[@]}"
	# do
	# ## Select fasta for each N subtype
	# # python ${scriptfolder}/select_by_meta.py --by "NSubtype=${nsub}" --fasta all.fasta --meta $metapath --output ${nsub}.fasta
	# ## Identity 90%, local alignment, 100% coverage for the shorter seq, any coverage for the longer seq, 20 threads
	# cd-hit -c 0.90 -G 0 -aS 1 -aL 0 -T 20 -i ${nsub}.fasta -o ${nsub}.09.cdhit.fasta
	# ## Adjust direction
	# mafft --retree 2 --maxiterate 2 --thread -1 --adjustdirection ${nsub}.09.cdhit.fasta >${alifolder}/${nsub}.fast.2it.msa 2>${alifolder}/${nsub}.fast.2it.msa.err
	# sed 's/-//g' ${alifolder}/${nsub}.fast.2it.msa >${alifolder}/${nsub}.fast.2it.temp
	# sed -i "/^$/d" ${alifolder}/${nsub}.fast.2it.temp
	# ## Cut gene regions from cds (for correct translation)
	# virulign refs/${nsub}.xml ${alifolder}/${nsub}.fast.2it.temp \
	# --exportKind GlobalAlignment --exportAlphabet Nucleotides --exportReferenceSequence yes --exportWithInsertions yes \
	# --progress yes --nt-debug ${alifolder}/Failed/${nsub} >${alifolder}/${nsub}.virualign 2>${alifolder}/${nsub}.virualign.err
	# done

## NOT USED
## Because a quarter of seqs do not have the id!
## Get coordinates from genbank??
# id was taken from meta.csv
# efetch -db nuccore  -format gbc -id AB049163 | xtract -pattern INSDSeq -group INSDFeature -if INSDFeature_key -equals gene -element INSDFeature_location
# cut -f1,60 /export/home/popova/workspace/birdflu/data/meta.csv >/export/home/popova/workspace/birdflu/data/NA_acc.list


## Merge protein fastas and remove duplicates 
# cat gisaid_epiflu_sequence.prot.fasta >all.prot.fasta.t
# cat gisaid_epiflu_sequence_since2015.prot.fasta >>all.prot.fasta.t
# faFilter -uniq all.prot.fasta.t all.prot.fasta


## tBLASTN proteins downloaded from gisaid to find the coordinates of the ORFs in nucleotide sequences %)
## Then translate these ORFs, align proteins with mafft and pal2nal
# for nsub in "${nsubs[@]}"
	# do
	## Select protein fasta for each N subtype
	# python ${scriptfolder}/select_by_meta.py --by "NSubtype=${nsub}" --fasta all.prot.fasta --meta $metapath --output ${nsub}.prot.fasta
	# ## Remove duplicates. Some ids have more than one seq. Nucleotide seqs can be different, but I assume that prots are the same (cRNA and DNA entries). 
	# # mv ${nsub}.1.cdhit.fasta ${nsub}.1.cdhit.fasta.t
	# # faFilter -uniq ${nsub}.1.cdhit.fasta.t ${nsub}.1.cdhit.fasta
	# # python ${scriptfolder}/find_coordinates.py --mode tblastn --alignment ${nsub}.prot.fasta \
	# # --fasta ${nsub}.1.cdhit.fasta --temp ${blastfolder}/temp \
	# # --output ${blastfolder}/${nsub}.blast
	## Instead of 2.cdhit.fasta, we should use fasta with adjusted direction
	## Remove gaps and remove duplicates:
	# sed 's/-//g' ${alifolder}/${nsub}.fast.2it.msa >${alifolder}/${nsub}.fast.2it.temp
	# sed -i "/^$/d" ${alifolder}/${nsub}.fast.2it.temp
	# faFilter -uniq ${alifolder}/${nsub}.fast.2it.temp ${alifolder}/${nsub}.fast.2it.uniq
	## tBlastn proteins against nucl fasta to get ORF coordinates:
	# python ${scriptfolder}/find_coordinates.py --mode tblastn --alignment ${nsub}.prot.fasta \
	# --fasta ${alifolder}/${nsub}.fast.2it.uniq --temp ${blastfolder}/temp \
	# --output ${blastfolder}/${nsub}.adj.blast
	## Get coding sequences which can be converted to correct proteins by transeq:
	# python ${scriptfolder}/cut_orfs.py --fasta ${alifolder}/${nsub}.fast.2it.uniq \
	# --tab ${blastfolder}/${nsub}.adj.blast --output ${alifolder}/${nsub}.coding
	# python ${scriptfolder}/transeq.py --input ${alifolder}/${nsub}.coding --output ${alifolder}/${nsub}.transeq.prot
	# mafft --anysymbol --thread -1 --maxiterate 10 --globalpair ${alifolder}/${nsub}.transeq.prot >${alifolder}/${nsub}.transeq.prot.slow10it
	# ## Convert protein alignment to nucleotide alignment
	# ## pal2nal does not understand Ns (but does understand dots)
	# perl -pi -e 's/N/./g unless $_ =~ "^>.*"' ${alifolder}/${nsub}.coding
	# ## pal2nal does not understand gapped codons (eg -TG), thus leading and trailing gaps in, respectively, leading and trailing codons are changed to dots
	# ## since we aligned seqs with reference proteins by codon-aware virulign, we do not expected any internal codons with gaps
	# python ${scriptfolder}/fix_gapfices.py --input ${alifolder}/${nsub}.coding --output ${alifolder}/${nsub}.coding.dotted
	# pal2nal.pl ${alifolder}/${nsub}.transeq.prot.slow10it ${alifolder}/${nsub}.coding.dotted -output fasta >${alifolder}/${nsub}.nucl.slow10it 2>${alifolder}/${nsub}.nucl.slow10it.err
	# done


## Merge intra-subtype alignments (without seed) 
>${alifolder}/concat.transeq.prot.slow10it
>${alifolder}/concat.coding.dotted
>${alifolder}/prot.subMSAtable
seqcount=1
for nsub in "${nsubs[@]}"
	do
		cat ${alifolder}/${nsub}.transeq.prot.slow10it >>${alifolder}/concat.transeq.prot.slow10it
		cat ${alifolder}/${nsub}.coding.dotted >>${alifolder}/concat.coding.dotted
		n=$(grep -c ">" ${alifolder}/${nsub}.transeq.prot.slow10it)
		t=$(($n+$seqcount))
		echo $(seq $seqcount $(($t-1)) | tr '\n' ' ' | sed 's/\s$//g') >>${alifolder}/prot.subMSAtable
		seqcount=$t
	done
## seems that --maxiterate does not work without --globalpair (or --localpair)
## without --anysymbol mafft changes * to a gap
## refinement is indispensable, since our sampling is very uneven (>50% seqs are from zaire ebolavirus) 
echo "Merging alignments for subtypes with mafft+pal2nal.."
mafft --anysymbol --thread -1 --maxiterate 10 --globalpair --merge ${alifolder}/prot.subMSAtable ${alifolder}/concat.transeq.prot.slow10it >${alifolder}/concat.transeq.prot.slow10it.slow10it
## We already prepared dotted file for pal2nal
pal2nal.pl ${alifolder}/concat.transeq.prot.slow10it.slow10it ${alifolder}/concat.coding.dotted -output fasta >${alifolder}/concat.transeq.nucl.slow10it.slow10it 2>${alifolder}/concat.transeq.nucl.slow10it.slow10it.err


## MANUALLY deleted and removed empty columns (jalview)
## N1.nucl.slow10it	EPI_ISL_273541	(gisaid link -> genbank:N2)
## N1.nucl.slow10it	EPI_ISL_290232	(N2)
## N1.nucl.slow10it	EPI_ISL_15724	(bad seq, several 1-nuc insertions and deletionss)s
## N2.nucl.slow10it	EPI_ISL_7035	("H6N2 was originally assigned by serotyping. After passage, sequence data indicate the genotype of the sample is H4N2,6.")
## N2.nucl.slow10it	EPI_ISL_167459	(closer to EPI_ISL_7035 than to any other N2 strain)
## N3.nucl.slow10it	EPI_ISL_267439	(H13N6 in genbank)
## N4.nucl.slow10it	EPI_ISL_9397	(H6N1)
## N4.nucl.slow10it	EPI_ISL_273544	(H3N6)
## N5.nucl.slow10it	EPI_ISL_24750	(blast -> genbank HM172178.1, H5N1)
## N5.nucl.slow10it	EPI_ISL_273611	(H7N3)
## N5.nucl.slow10it	EPI_ISL_273549	(H10N7)
## N5.nucl.slow10it	EPI_ISL_273598	(H4N6 Two seqs, second one is H10N5 but we have the first)
## N6.nucl.slow10it	EPI_ISL_7036	(mixed H6N1,6)
## N8.nucl.slow10it	EPI_ISL_7034
## N8.nucl.slow10it	EPI_ISL_257145
## N8.nucl.slow10it	EPI_ISL_75721
## N8.nucl.slow10it	EPI_ISL_267354
## N8.nucl.slow10it	EPI_ISL_273594
## N9.nucl.slow10it	EPI_ISL_4057
## N9.nucl.slow10it	EPI_ISL_267228


