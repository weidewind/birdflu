from Bio import SeqIO
from Bio.Seq import Seq
import optparse

parser = optparse.OptionParser()
parser.add_option('--fasta', help='', type='str')
parser.add_option('--start', help='', type='int')
parser.add_option('--stop', help='', type='int')
parser.add_option('--output', help='', type='str')

options, args = parser.parse_args()

fasta = SeqIO.parse(options.fasta, "fasta")
for r in fasta:
	with open(options.output, "w") as out:
		out.write(">" + r.id + " " + str(options.start) + ":" + str(options.stop) + "\n")
		out.write(str(r.seq[options.start-1:options.stop]) + "\n")

