
import os
import sys
from collections import defaultdict
from Bio import SeqIO
import argparse

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Extract sequence for segments from assembly.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--scaffolds_fasta', help='Multi-FASTA of all scaffolds in analysis.', required=True)
	parser.add_argument('-s', '--segments_listing', help='Concatenated Segment_Results.txt files from delineateSegmentsOnReference.py.', required=True)
	args = parser.parse_args()
	return args

myargs = create_parser()
scaffolds_fasta = myargs.scaffolds_fasta
segments_listing = myargs.segments_listing

scaff_seqs = {}
scaff_lens = {}
with open(scaffolds_fasta) as oaf:
	for rec in SeqIO.parse(oaf, 'fasta'):
		scaff_seqs[rec.id] = str(rec.seq)
		scaff_lens[rec.id] = len(str(rec.seq))

with open(segments_listing) as odf:
	for line in odf:
		line = line.strip()
		ref_scaff, start, end, scaffolds = line.split('\t')
		scaffolds = scaffolds.split(', ')
		el_id = ref_scaff + '|' + start + '|' + end
		subseq = ""
		seq = scaff_seqs[ref_scaff]
		reflen = scaff_lens[ref_scaff]
		start = int(start)
		end = int(end)
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

		if subseq != None:
			print('>%s|%d|%d\n%s' % (ref_scaff, start, end, subseq))
		else:
			sys.stderr.write("Issue extracting sequence for %s...\n" % ref_scaff)