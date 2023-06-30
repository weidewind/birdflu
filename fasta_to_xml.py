import optparse
from Bio import SeqIO
import re

parser = optparse.OptionParser()
parser.add_option('--input', help='', type='str')
parser.add_option('--output', help='msa file', type='str')


options, args = parser.parse_args()
fasta = SeqIO.parse(options.input, "fasta")
for record in fasta:
	res = re.search(r'([0-9]+)\.\.([0-9]+)', record.description)
	print(record.id)
	print(res)
	with open(options.output + record.id + ".xml", "w") as out:
		out.write('<?xml version="1.0" encoding="UTF-8"?>\n')
		out.write('<orf name="' + record.id + '" referenceSequence="' + str(record.seq) + '" >\n')
		## res[0] is the whole captured regex
		out.write('<protein abbreviation="' + record.id + '" startPosition="' + res[1] + '" stopPosition="'+ res[2] +'" />\n')
		out.write('</orf>')