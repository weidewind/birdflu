## For protein-coding nucleotide alignments
## Translates -ATATGATGAT- as NMMN

from Bio import SeqIO
from Bio.Seq import Seq
import optparse

parser = optparse.OptionParser()
parser.add_option('--input', help='', type='str')
parser.add_option('--N', help='', type='str', default='N')
parser.add_option('--output', help='', type='str')

options, args = parser.parse_args()

def translate(dna):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        '---':'-'
        }
    ltable = {k.lower():v for k,v in table.items()}
    table.update(ltable)
    protein = []
    end = len(dna) - (len(dna) %3) - 1
    for i in range(0,end,3):
        codon = dna[i:i+3]
        if codon in table:
            aa = table[codon]
            protein.append(aa)
        else:
            protein.append("X")
    return ("".join(protein))


records = SeqIO.parse(options.input, "fasta")
with open(options.output, "w") as out:
	for r in records:
		out.write(">" + r.id + "\n")
		out.write(translate(r.seq) + "\n")


