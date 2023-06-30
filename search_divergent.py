from scipy import stats
import optparse
import pandas as pd
from Bio import SeqIO, AlignIO
from collections import deque, Counter
import numpy as np
import sys
import divergent_table as dtable
 
parser = optparse.OptionParser()
parser.add_option('--table', help='path', type='str')
parser.add_option('--window', help='', type='int')
parser.add_option('--in_cons_thr', help='frequency of the consensus symbol in my group must be at least in_cons_thr', type='float')
parser.add_option('--out_cons_thr', help='frequency of the cconsensus symbol (taken from my group) outside of my group must be no more than out_cons_thr', type='float')
parser.add_option('--output', help='', type='str')
parser.add_option('--groups', help='comma-delimited list of subtypes to compare to', type='str')

options, args = parser.parse_args()

table = pd.read_csv(options.table, sep = "\t", comment='#')

length = table.shape[0]
deck = deque()

# default_groups = ["N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9"]
# if options.groups:
# 	groups = options.groups.split(",")
# 	for g in groups:
# 		if g not in default_groups:
# 			sys.exit('No such N subtype:' + g + "! Please provide comma-delimited list of neuraminidase subtypes\n")
# else:
# 	groups = default_groups

# outgroup_consfreq_colnames = [g + "_cons_freq" for g in groups]


consfreq_colnames = dtable.get_cfreq_colnames(table.columns)
if options.groups:
	outgroup_consfreq_colnames = [g + "_" + dtable.cfreq_suffix() for g in options.groups.split(",")]
	for og in outgroup_consfreq_colnames:
		if og not in consfreq_colnames:
			sys.exit('No such column:' + og + " in divergence table! Please provide comma-delimited list of neuraminidase subtypes\n")
else:
	outgroup_consfreq_colnames = consfreq_colnames

print(outgroup_consfreq_colnames)

cfreq_suffix = dtable.cfreq_suffix()
gfreq_suffix = dtable.gfreq_suffix()

mg_deq = deque()
og_deq = deque()
mg_counter = 0
og_counter = 0
mg_list = []
og_list = []

og_count_deq = deque()
og_count_counter = 0
og_count_list = []

gap_deq = deque()
gap_counter = 0
gap_list = []

## number of positions in window, where Freq_consensus_symbol in my group >= in_cons_thr
def mg_count(row, mg_deq, mg_counter):
	if row[cfreq_suffix] >= options.in_cons_thr:
		mg_deq.append(1)
		mg_counter += 1
	else:
		mg_deq.append(0)
	if len(mg_deq) > options.window:
		mg_counter -= mg_deq.popleft()
	return mg_deq, mg_counter

## number of positions in window, where Freq_consensus_symbol (taken from my group) in each of the other groups <= out_cons_thr
def og_count(row, og_deq, og_counter):
	if all([row[gcf] <= options.out_cons_thr for gcf in outgroup_consfreq_colnames]):
		og_deq.append(1)
		og_counter += 1
	else:
		og_deq.append(0)
	if len(og_deq) > options.window:
		og_counter -= og_deq.popleft()
	return og_deq, og_counter

## total Freq_consensus_symbol in my group
def mg_total(row, mg_deq, mg_counter):
	mg_deq.append(row[cfreq_suffix])
	mg_counter += row[cfreq_suffix]
	if len(mg_deq) > options.window:
		mg_counter -= mg_deq.popleft()
	return mg_deq, mg_counter

## total Freq_consensus_symbol (taken from my group) in the other groups
def og_total(row, og_deq, og_counter):
	for gcf in outgroup_consfreq_colnames:
		og_deq.append(row[gcf])
		og_counter += row[gcf]
	if len(og_deq) > options.window:
		og_counter -= og_deq.popleft()
	return og_deq, og_counter

## total Freq_gaps in my group
def gap_total(row, gap_deq, gap_counter):
	gap_deq.append(row[gfreq_suffix])
	gap_counter += row[gfreq_suffix]
	if len(gap_deq) > options.window:
		gap_counter -= gap_deq.popleft()
	return gap_deq, gap_counter


## TODO exclude windows with too many gaps (missing indices)!
## Now we just ignore missing indices
for index, row in table.iterrows():
	
	# mg_deq, mg_counter = mg_count(row, mg_deq, mg_counter)
	# og_deq, og_counter = og_count(row, og_deq, og_counter)

	mg_deq, mg_counter = mg_total(row, mg_deq, mg_counter)
	og_deq, og_counter = og_total(row, og_deq, og_counter)
	og_count_deq, og_count_counter = og_count(row, og_count_deq, og_count_counter)
	gap_deq, gap_counter = gap_total(row, gap_deq, gap_counter)

	mg_list.append(mg_counter/len(mg_deq))
	og_list.append(og_counter/len(og_deq))
	og_count_list.append(og_count_counter/len(og_count_deq))
	gap_list.append(gap_counter/len(gap_deq))


with open(options.output, "w") as out:
	out.write("\t".join(["start_ind", "end_ind", "GapFreq", "Nc_mygroup", "Nnc_outgroups", "FreqFreq<" + str(options.out_cons_thr)]) + "\n")
	for index, row in table.iterrows():
		start_row_ind = max(0, index - options.window + 1)
		out.write("\t".join([str(table.iloc[start_row_ind]["ind"]), str(row["ind"]), str(gap_list[index]), str(mg_list[index]), str(og_list[index]), str(og_count_list[index])]) + "\n")

	## print options
	for opt,val in vars(options).items():
		out.write("## " + opt + "\t" + "'" + str(val) + "'" +"\n")
