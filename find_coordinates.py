from Bio import SeqIO
from Bio.Seq import Seq
import optparse
import subprocess
import os
from multiprocessing import Manager
from multiprocessing.pool import Pool

parser = optparse.OptionParser()
parser.add_option('--alignment', help='', type='str')
parser.add_option('--fasta', help='', type='str')
parser.add_option('--temp', help='path to folder', type='str')
parser.add_option('--output', help='', type='str')
parser.add_option('--mode', help='', type='str', default='blastn')

options, args = parser.parse_args()


def blastn(q, fmt, lock):
	if q.id in subjects:
		query_path = os.path.join(options.temp, q.id + ".query.fasta")
		with open(query_path, "w") as out:
			out.write(">" + q.id + "\n" + str(q.seq))
		subj_path = os.path.join(options.temp, q.id + ".subject.fasta")
		with open(subj_path, "w") as out:
			out.write(">" + q.id + "\n" + str(subjects[q.id].seq))
		blout_path = os.path.join(options.temp, q.id + ".blast.out")
		subprocess.run([options.mode, '-query', query_path, '-subject', subj_path, '-outfmt', fmt, '-out', blout_path])
		with lock:
			with open(options.output, 'a+') as blastout:
				add = open(blout_path, 'r')
				blastout.write(add.read())



# protect the entry point
if __name__ == '__main__':
	queries = SeqIO.parse(options.alignment, "fasta")
	subjects = SeqIO.to_dict(SeqIO.parse(options.fasta, "fasta"))
	#	subprocess.run(['cp', '/export/home/popova/workspace/filoviridae/output/blast/outfmt6.header', options.output])
	vals = ["qseqid", "sseqid", "pident", "qlen", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
	fmt = "6 " + " ".join(vals)
	with open(options.output, "w") as out:
		out.write("\t".join(vals) + "\n")
	with Manager() as manager:
		lock = manager.Lock()
		with Pool(10) as pool:
			items = [(q, fmt, lock) for q in queries]
			result = pool.starmap_async(blastn, items)
			result.wait()


