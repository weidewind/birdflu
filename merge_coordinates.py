import pandas as pd
import optparse
import os

parser = optparse.OptionParser()
parser.add_option('--blast_folder', help='', type='str')
parser.add_option('--output', help='', type='str')
parser.add_option('--genes', help='comma-separated list', type='str', default='NP,VP35,VP40,GP,VP30,VP24,L')
parser.add_option('--ext', help='eg, .adj.blast', type='str')

options, args = parser.parse_args()

# genes = ["NP", "VP35", "VP40", "GP", "VP30", "VP24", "L"]
# colnames = ["id"]
# colnames.extend([gene + "_start" for gene in genes])
# colnames.extend([gene + "_end" for gene in genes])
# table = pd.DataFrame(columns = colnames)
# table.set_index('id', inplace=True)

def gene_df(gene):
	print(gene)
	gtable = pd.read_csv(os.path.join(options.blast_folder, gene + options.ext), sep="\t")
	print(gtable[0:5])
	#gtable = gtable[(gtable['qstart'] == 1) & (gtable['qend'] == gtable['length'])][['qseqid', 'sstart', 'send']]
	full = gtable[(gtable['qlen'] == gtable['length'])]  
	short = gtable[~gtable.qseqid.isin(full.qseqid.tolist())]
	short = short.sort_values(['length']).drop_duplicates(['qseqid'], keep='last')
	short['sstart'] = short['sstart'] - short['qstart'] + 1
	short['send'] = short['send'] + short['qlen'] - short['qend']
	gtable = full.append(short)[['qseqid', 'sstart', 'send']]
	gtable.rename(columns={'qseqid':'id', 'sstart':gene + "_start", 'send':gene + "_end"}, inplace=True)
	gtable.set_index('id', inplace=True, verify_integrity=True)
	return(gtable)

genes = options.genes.split(",")
table = gene_df(genes.pop(0))
for gene in genes:
	gtable = gene_df(gene)
	table[gene + "_start"] = gtable[gene + "_start"]
	table[gene + "_end"] = gtable[gene + "_end"]
	#table = table.combine_first(gtable)
	#table = pd.merge(table, gtable, on = 'id', how = 'outer')

table.to_csv(options.output, index = True, sep = "\t")
