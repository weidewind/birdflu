

datafolder=/export/home/popova/workspace/birdflu/data
alifolder=/export/home/popova/workspace/birdflu/output/aligned_2
scriptfolder=/export/home/popova/workspace/birdflu/scripts
metapath=${datafolder}/meta_short.csv
cd $datafolder
nsubs=(N9) #N4 N5 N1 N2 N3 N6 N7 N8 

# ## Make meta_short.csv
# python xls_to_csv.py --xls gisaid_epiflu_isolates.xls --csv gisaid_epiflu_isolates.csv
# python xls_to_csv.py --xls gisaid_epiflu_isolates_since2015.xls --csv gisaid_epiflu_isolates_since2015.csv
# cat gisaid_epiflu_isolates.csv >meta.csv
# tail -n +2 gisaid_epiflu_isolates_since2015.csv >>meta.csv
# cut -f1,13,14,15,16,17,18,19,26 meta.csv >${metapath}
# python add_Nsubtype.py

# ## Merge fastas
# cat gisaid_epiflu_sequence_onlyid.fasta >all.fasta
# cat gisaid_epiflu_sequence_onlyid_since2015.fasta >>all.fasta
## MANUALLY delete 
## >EPI_ISL_11050875
## achickenbangladeshcatgaatcc
# sed -i 's/achickenbangladesh//g' ../data/all.fasta

# ## Create xmls for VIRULIGN from the reference file
# ##  References were downloaded from genbank quite randomly (avian, complete cds)
# python fasta_to_xml.py --input refs/refs.nucl.fasta --output refs/


# ## Count entries attributed to different subtypes.
# ## Not used, I just wanted to look at them
# cut -f3 meta_short.csv | sort | uniq -c >Nsubtypes.list
# cut -f2 meta_short.csv | sort | uniq -c >subtypes.list



## Virulign needs path export: screen -U -S virnp env LD_LIBRARY_PATH=$LD_LIBRARY_PATH 
for nsub in "${nsubs[@]}"
	do
	# ## Select fasta for each N subtype
	# python ${scriptfolder}/select_by_meta.py --by "NSubtype=${nsub}" --fasta all.fasta --meta $metapath --output ${nsub}.fasta
	# ## Identity 100%, local alignment, 100% coverage for the shorter seq, any coverage for the longer seq, 20 threads
	# cd-hit -c 1 -G 0 -aS 1 -aL 0 -T 20 -i ${nsub}.fasta -o ${nsub}.1.cdhit.fasta
	## Adjust direction
	# mafft --retree 2 --maxiterate 2 --thread -1 --adjustdirection ${nsub}.1.cdhit.fasta >${alifolder}/${nsub}.fast.2it.msa 2>${alifolder}/${nsub}.fast.2it.msa.err
	## Cut gene regions from cds (for correct translation)
	virulign refs/${nsub}.xml ${alifolder}/${nsub}.fast.2it.msa \
	--exportKind GlobalAlignment --exportAlphabet Nucleotides --exportReferenceSequence yes --exportWithInsertions yes \
	--progress yes --nt-debug ${alifolder}/Failed/${nsub} >${alifolder}/${nsub}.virualign 2>${alifolder}/${nsub}.virualign.err
	## Align proteins
	python ${scriptfolder}/transeq.py --input ${alifolder}/${nsub}.virualign --output ${alifolder}/${nsub}.virualign.prot
	mafft --anysymbol --thread -1 --maxiterate 10 --globalpair ${alifolder}/${nsub}.virualign.prot >${alifolder}/${nsub}.virualign.prot.slow10it
	## Convert protein alignment to nucleotide alignment
	## pal2nal does not understand Ns (but does understand dots)
	perl -pi -e 's/N/./g unless $_ =~ "^>.*"' ${alifolder}/${nsub}.virualign
	## pal2nal does not understand gapped codons (eg -TG), thus leading and trailing gaps in, respectively, leading and trailing codons are changed to dots
	## since we aligned seqs with reference proteins by codon-aware virulign, we do not expected any internal codons with gaps
	python ${scriptfolder}/fix_gapfices.py --input ${alifolder}/${nsub}.virualign --output ${alifolder}/${nsub}.virualign.dotted
	pal2nal.pl ${alifolder}/${nsub}.virualign.prot.slow10it ${alifolder}/${nsub}.virualign.dotted -output fasta >${alifolder}/${nsub}.virualign.nucl.slow10it 2>${alifolder}/${nsub}.virualign.nucl.slow10it.err
	done
