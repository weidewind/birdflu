from scipy import stats
import optparse
import pandas as pd
from Bio import SeqIO, AlignIO
from collections import deque, Counter
import numpy as np
import sys
import divergent_table as dtable
 
parser = optparse.OptionParser()
parser.add_option('--msa', help='reference fasta', type='str')
parser.add_option('--meta', help='', type='str')
parser.add_option('--pivot', help='subtype to compare the others to', type='str')
parser.add_option('--numeration_id', help='', type='str')
parser.add_option('--output', help='', type='str')

options, args = parser.parse_args()

metafile = options.meta
meta = pd.read_csv(metafile, sep = "\t")
group_dict =  dict(zip(meta["Isolate_Id"], meta["NSubtype"]))

## read() assumes there is only one msa in the file. parse() returns a list of msas! 
msa = AlignIO.read(options.msa, "fasta")
groups = ["N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9"]
for r in msa:
	r.id = r.id.split("/")[0]
	## If numeration follows options.numeration_id, we need numeration_seq
	if options.numeration_id and r.id == options.numeration_id:
		numeration_seq = r.seq
if numeration_seq is None:
	sys.exit('No ' + options.numeration_id + " found in msa " + options.msa + "!\n")

msa_dict = {}
## k: group
## v: [msa records assigned to the group]
width_dict = {}
## k: group
## v: number of seqs
for g in groups:
	seqs = [r for r in msa if group_dict[r.id] == g]
	msa_dict[g] = seqs
	width_dict[g] = len(seqs)


my_group = options.pivot
other_groups = [g for g in groups if g != my_group]


## For every msa column, print consensus char and freqs of this char in all groups
counter_dict = {ind:{} for ind in range(msa.get_alignment_length())}
## k: index
## v: dict{group: dict{symbol:count}}
with open(options.output, "w") as out:

	## Header
	# out.write("\t".join(["ind", "cons", "cons_freq", "gap_freq"]))
	# for g in other_groups:
	# 	out.write("\t" + g + "_cons_freq" + "\t" + g + "_gap_freq")
	# out.write("\n")
	out.write(dtable.get_header(other_groups) + "\n")

	seq_gaps = 0
	## seq_gaps counts gaps in numeration_seq
	for ind in range(msa.get_alignment_length()):
		if options.numeration_id and numeration_seq[ind] == "-":
			seq_gaps += 1
		## Consensus symbol: if several symbols have the same occurence, just take the first
		## REMOVED If the most frequent symbol is gap, take the next most frequent symbol
		## REMOVED If gap is the only symbol, skip the column
		mg_col = [r.seq[ind] for r in msa_dict[my_group]]
		mg_symbol_count = Counter(mg_col)
		counter_dict[ind][my_group] = mg_symbol_count
		mg_gap_count = mg_symbol_count.get("-", 0)
		# if mg_gap_count == width_dict[my_group]:
		# 	continue
		mg_cons,mg_cons_count = mg_symbol_count.most_common(1)[0]
		# if mg_cons == "-":
		# 	mg_cons,mg_cons_count = mg_symbol_count.most_common(2)[1]
		# mg_cons_freq = mg_cons_count/(width_dict[my_group] - mg_gap_count)
		mg_cons_freq = mg_cons_count/width_dict[my_group]
		mg_gap_freq = mg_gap_count/width_dict[my_group]

		if options.numeration_id:
			index = ind + 1 - seq_gaps
		else:
			index = 0
		out.write("\t".join([str(n) for n in [index, mg_cons, mg_cons_freq, mg_gap_freq]]))
		for g in other_groups:
			col = [r.seq[ind] for r in msa_dict[g]]
			symbol_count = Counter(col)
			counter_dict[ind][g] = symbol_count
			gap_count = symbol_count.get("-", 0)
			# if gap_count == width_dict[g]:
			# 	cons_freq = 0
			# else:
			# 	cons_freq = symbol_count[mg_cons]/(width_dict[g] - gap_count)
			cons_freq = symbol_count[mg_cons]/width_dict[g]
			gap_freq = gap_count/width_dict[g]
			out.write("\t" + "\t".join([str(n) for n in [cons_freq, gap_freq]]))
		out.write("\n")

	for opt,val in vars(options).items():
		out.write("## " + opt + "\t" + str(val) + "\n")