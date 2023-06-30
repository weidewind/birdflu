

datafolder=/export/home/popova/workspace/birdflu/data
outfolder=/export/home/popova/workspace/birdflu/output
scriptfolder=/export/home/popova/workspace/birdflu/scripts
metapath=${datafolder}/meta_short.csv
cd $datafolder
nsubs=(N4 N5 N1 N2 N3 N6 N7 N8 N9)

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
## TODO remove duplicates

faFilter -minSize=1300 -maxN=70 all.fasta all.SizeN.filtered.fasta
grep ">" all.fasta>all.fasta.ids
grep ">" all.SizeN.filtered.fasta >all.SizeN.filtered.fasta.ids
sort all.fasta.ids >all.fasta.ids.sorted.t
tr -d '\015' <all.fasta.ids.sorted.t >all.fasta.ids.sorted
sort all.SizeN.filtered.fasta.ids >all.SizeN.filtered.fasta.ids.sorted
comm -23 all.fasta.ids.sorted all.SizeN.filtered.fasta.ids.sorted | tr -d '>' >removed.SizeN.list
cat removed.SizeN.list >remove.list
cut -f1 removeMistyped.list >>remove.list
faFilter -namePatList=remove.list -v all.fasta all.filtered.fasta

## Filter fasta and make a remove.idlist file with ids that need to be removed.
## We can use it in downstream analysis


# ## Create xmls for VIRULIGN from the reference file
# ##  References were downloaded from genbank quite randomly (avian, complete cds)
# python ${scriptfolder}/fasta_to_xml.py --input refs/refs.nucl.fasta --output refs/


# ## Count entries attributed to different subtypes.
# ## Not used, I just wanted to look at them
# cut -f3 meta_short.csv | sort | uniq -c >Nsubtypes.list
# cut -f2 meta_short.csv | sort | uniq -c >subtypes.list


## Select fasta for each N subtype
## Virulign needs path export: screen -U -S virnp env LD_LIBRARY_PATH=$LD_LIBRARY_PATH 
for nsub in "${nsubs[@]}"
	do
	# python ${scriptfolder}/select_by_meta.py --by "NSubtype=${nsub}" --fasta all.fasta --meta $metapath --output ${nsub}.fasta
	# ## Identity 100%, local alignment, 100% coverage for the shorter seq, any coverage for the longer seq, 20 threads
	# cd-hit -c 1 -G 0 -aS 1 -aL 0 -T 20 -i ${nsub}.fasta -o ${nsub}.1.cdhit.fasta
	mafft --retree 2 --maxiterate 100 --thread -1 --adjustdirection ${nsub}.1.cdhit.fasta >${outfolder}/aligned/${nsub}.fast.100it.msa 2>${outfolder}/aligned/${nsub}.fast.100it.msa.err
	virulign refs/${nsub}.xml ${outfolder}/aligned/${nsub}.fast.100it.msa \
	--exportKind GlobalAlignment --exportAlphabet Nucleotides --exportReferenceSequence yes --exportWithInsertions yes \
	--progress yes --nt-debug ${outfolder}/aligned/Failed/${nsub} >${outfolder}/aligned/${nsub}.virualign 2>${outfolder}/aligned/${nsub}.virualign.err
	done
