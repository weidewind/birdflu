import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import optparse
import os

parser = optparse.OptionParser()
parser.add_option('--fasta', help='', type='str')
parser.add_option('--tab', help='', type='str')
parser.add_option('--output', help='', type='str')

options, args = parser.parse_args()

tab = pd.read_csv(options.tab, sep = "\t")
print(tab[0:5])
fasta_dict = SeqIO.to_dict(SeqIO.parse(options.fasta, "fasta"))

with open(options.output, "w") as out:
	for index, row in tab.iterrows():
		if row["qlen"] == row["length"]:
			out.write(">" + row["qseqid"] + "\n")
			seq = str(fasta_dict[row["qseqid"]].seq)
			out.write(seq[(row["sstart"]-1):row["send"]] + "\n")
