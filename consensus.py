import optparse
from Bio import SeqIO, AlignIO, Align
from Bio.Align import AlignInfo
 
parser = optparse.OptionParser()
parser.add_option('--msa', help='reference fasta', type='str')
parser.add_option('--name', help='', type='str')
parser.add_option('--output', help='', type='str')

options, args = parser.parse_args()


## read() assumes there is only one msa in the file. parse() returns a list of msas! 
msa = AlignIO.read(options.msa, "fasta")
summary = AlignInfo.SummaryInfo(msa)
with open(options.output, "w") as out:
	out.write(">" + options.name + "\n" + str(summary.gap_consensus(threshold=0.5, ambiguous='X')) + "\n")
