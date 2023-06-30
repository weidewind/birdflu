## if a seq has gap prefix, change 0-2 gaps preceding the seq to . (dot)
## so that all codons in this ORF consisted only of letters and dots
## eg, ----TG -> ---.TG 
## Same for gap suffix
## eg, --TGACG----- -> ..TGACG..---
## (needed for pal2nal)

from Bio import SeqIO
from Bio.Seq import Seq
import re
import optparse

parser = optparse.OptionParser()
parser.add_option('--input', help='', type='str')
parser.add_option('--output', help='', type='str')

options, args = parser.parse_args()

records = SeqIO.parse(options.input, "fasta")
with open(options.output, "w") as out:
	for r in records:
		print(r.id)
		gappref = re.search(r'^-+', str(r.seq))
		new_gappref = ""
		dotpref  = ""
		if gappref:
			gappref = gappref.group()
			res = len(gappref)%3
			new_gappref = gappref
			if res > 0:
				new_gappref = "---" * int(len(gappref)/3)
				dotpref = "." * res
		gapsuf = re.search(r'-+$', str(r.seq))
		new_gapsuf = ""
		dotsuf = ""
		if gapsuf:
			gapsuf = gapsuf.group()
			head = r.seq[0:(len(r.seq)-len(gapsuf))] 
			print(head)
			res = len(head)%3
			if res > 0:
				dotsuf = "." * (3 - res)
				new_gapsuf = "-" * (len(gapsuf) - (3 -res))
		new_seq = new_gappref + dotpref + str(r.seq[(len(new_gappref) + len(dotpref)):(len(r.seq) -len(dotsuf) - len(new_gapsuf))]) + dotsuf + new_gapsuf
		out.write(">" + r.id + "\n" + new_seq + "\n")
