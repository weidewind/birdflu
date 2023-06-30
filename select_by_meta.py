import optparse
import pandas as pd
from Bio import SeqIO
import re

parser = optparse.OptionParser()
parser.add_option('--fasta', help='', type='str')
parser.add_option('--meta', help='', type='str')
parser.add_option('--by', help='eg, protein=nucleoprotein,Subtype=A / H9N2', type='str')
parser.add_option('--output', help='', type='str')

options, args = parser.parse_args()

metafile = options.meta
meta = pd.read_csv(metafile, sep = "\t")
features = {f.split("=")[0]:f.split("=")[1] for f in options.by.split(",")}
print(features)

## https://stackoverflow.com/questions/35882501/selecting-data-from-pandas-dataframe-based-on-criteria-stored-in-a-dict
mask = pd.DataFrame([meta[f] == val for f, val in features.items()]).T.all(axis=1)
ids = set(meta[mask]["Isolate_Id"].tolist())

records = SeqIO.parse(options.fasta, "fasta")
with open(options.output, "w") as out:
	for r in records:
		if r.id in ids:
			out.write(">" + r.id + "\n")
			out.write(str(r.seq) + "\n")

