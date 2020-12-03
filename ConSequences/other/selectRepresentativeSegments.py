import os
import sys
from collections import defaultdict
import argparse

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Parse CD-HIT results and select representatives from clusters of analogous segments.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-s', '--segments_listing', help='Concatenated Segment_Results.txt files from delineateSegmentsOnReference.py.', required=True)
	parser.add_argument('-c', '--cdhit_clustering', help='CD-HIT based clustering.', required=True)
	args = parser.parse_args()
	return args

myargs = create_parser()
segments_listing = myargs.segments_listing
cdhit_clustering = myargs.cdhit_clustering

# Note, for Pironti and Salamzade et al. 2020, we selected representative segments based on maximizing the number
# of samples; however, here we implement based on maximizing the number of scaffolds.

core_scaffolds = defaultdict(set)
de_info = {}
with open(segments_listing) as odf:
	for line in odf:
		line = line.strip()
		ref_scaff, start, stop, scaffolds = line.split('\t')
		scaffolds = scaffolds.split(', ')
		el_id = ref_scaff + '|' + start + '|' + stop
		core_scaffolds[el_id] = scaffolds
		de_info[el_id] = line

cluster_mems = defaultdict(set)
cluster = None
with open(cdhit_clustering) as ocf:
	for line in ocf:
		line = line.strip()
		if line.startswith(">"):
			cluster = line[1:]
		else:
			mem = line.split()[2][1:-3]
			memlen = int(line.split()[1][:-3])
			cluster_mems[cluster].add(tuple([mem, memlen]))

for c in cluster_mems:
	mems = cluster_mems[c]

	max_core_size = [[], 0]
	for m in mems:
		if len(core_scaffolds[m[0]]) > max_core_size[1]: max_core_size = [[m], len(core_scaffolds[m[0]])]
		elif len(core_scaffolds[m[0]]) == max_core_size[1]: max_core_size[0].append(m)

	best_m = [0, 0]
	for m in sorted(max_core_size[0]):
		if m[1] > best_m[1]:
			best_m = [m[0], m[1]]

	print(de_info[best_m[0]])