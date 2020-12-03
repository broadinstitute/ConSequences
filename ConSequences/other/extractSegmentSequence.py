import os
import sys
from collections import defaultdict
from Bio import SeqIO
import argparse

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Extract sequence for segment from assembly.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--scaffolds_fasta', help='Multi-FASTA of all scaffolds in analysis.', required=True)
	parser.add_argument('-r', '--scaffold', help='Scaffold identifier.', required=True)
	parser.add_argument('-s', '--start', help='Start position.', required=True)
	parser.add_argument('-e', '--end', help='End position.', required=True)
	args = parser.parse_args()
	return args

myargs = create_parser()
scaffolds_fasta = myargs.scaffolds_fasta
scaffold = myargs.scaffold
start = myargs.start
end = myargs.end

subseq = ""
with open(scaffolds_fasta) as oaf:
	for rec in SeqIO.parse(oaf, 'fasta'):
		if rec.id == scaffold:
			seq = str(rec.seq)
			if start > end:
				pos_range = set(range(start-1, reflen))
				for p, bp in enumerate(seq):
					if p in pos_range: subseq += bp
				pos_range = set(range(0, end))
				for p, bp in enumerate(seq):
					if p in pos_range: subseq += bp

			else:
				pos_range = set(range(start-1, end))
				for p, bp in enumerate(seq):
					if p in pos_range: subseq += bp

if subseq != "":
	print('>%s|%d|%d\n%s' % (scaffold, start, end, subseq))
else:
	sys.stderr.write("Issue extracting sequence! Exiting now...\n")
	sys.exit(1)