#!/usr/bin/env python

### Program: delineateSegmentsOnReference.py
### Author: Rauf Salamzade
### The Broad Institute of MIT and Harvard
### Earl Lab / Bacterial Genomics Group

"""
Special thanks to: Abigail Manson on feedback and improvement suggestions around workflow and Colin Worby for
solution to append sliding-window's worth of sequence from end to beginning of scaffold to account for plasmid
circularity/completeness.
"""

import os
import sys
from time import sleep
import argparse
from ConSequences import ConSequences

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: delineateSegmentsOnReference.py
	Author: Rauf Salamzade
	The Broad Institute of MIT and Harvard
	Earl Lab / Bacterial Genomics Group

	This program will take a reference scaffold and identify conserved and contiguous segments it shares with one or more
	query scaffolds. If facing difficulties, please raise issues on the github page: https://github.com/broadinstitute/consequences
	""", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('-r', '--ref_fasta', type=str, help="FASTA (single entry) for reference scaffold upon which to call windows on.", required=True)
	parser.add_argument('-q', '--query_multifasta', type=str, help="Multi-FASTA for query scaffolds to use in search (should include reference scaffold as well).", required=True)
	parser.add_argument('-o', '--outdir', type=str, help="", required=True)
	parser.add_argument('-rc', '--ref_circular', action='store_true', help="Is the reference scaffold circular/complete?", default=False, required=False)
	parser.add_argument('-qc', '--query_circular', type=str, help="A file specifying which of the query scaffolds are circular/complete. Each query identifier should be a separate line and match an identifier in the query Multi-FASTA file", default="", required=False)
	parser.add_argument('-w', '--scan_window_size', type=int, help='length of windows to use.', required=False, default=10000)
	parser.add_argument('-s', '--scan_slide_step', type=int, help='granularity of sliding.', required=False, default=100)
	parser.add_argument('-c', '--matches_percentage', type=float, help='cutoff for number of matches in common.', required=False, default=99.0)
	parser.add_argument('-f', '--unfilter_segs', action='store_true', help="Do not filter delineated segments for those representative of the sample set of a specific window.", required=False, default=False)
	args = parser.parse_args()
	return args

def main():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE REQUIRED INPUTS
	"""

	myargs = create_parser()

	ref_fasta_file = os.path.abspath(myargs.ref_fasta)
	query_multifasta_file = os.path.abspath(myargs.query_multifasta)
	outdir = os.path.abspath(myargs.outdir) + '/'

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n ")
		sleep(5)
	else:
		os.system('mkdir %s' % outdir)

	try:
		assert(ConSequences.is_fasta(ref_fasta_file) and ConSequences.is_fasta(query_multifasta_file))
	except:
		raise RuntimeError("Input FASTA files are invalid. Raising error.")

	"""
	PARSE OPTIONAL INPUTS
	"""

	ref_circular = myargs.ref_circular
	query_circular_file = myargs.query_circular
	scan_window_size = myargs.scan_window_size
	scan_slide_step = myargs.scan_slide_step
	matches_percentage = myargs.matches_percentage
	unfilter_segs = myargs.unfilter_segs

	query_circular = set([])
	if os.path.isfile(query_circular_file):
		query_circular_file = os.path.abspath(query_circular_file)
		with open(query_circular_file) as oqcf:
			for line in oqcf:
				query_circular.add(line.strip())

	"""
	START WORKFLOW
	"""

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = ConSequences.createLoggerObject(log_file)

	# Step 0: Log input arguments and update reference and query FASTA files.

	logObject.info("Saving parameters for future provedance.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [ref_fasta_file, query_multifasta_file, outdir, ref_circular, query_circular_file,
					   scan_slide_step, scan_window_size, matches_percentage, unfilter_segs]
	parameter_names = ["Reference FASTA file", "Query Multi-FASTA file", "Output Directory", "Reference is Circular?",
						"Query Scaffolds which are Circular file", "Scan Window Size", "Scan Window Step",
					    "Matches Percentage", "Unfilter Segments"]
	ConSequences.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	logObject.info("Running preparation of FASTA inputs to account for circularity/completeness of scaffolds.")
	ref_fasta_file, query_multifasta_file = ConSequences.run_prep(ref_fasta_file, query_multifasta_file, scan_window_size, ref_circular, query_circular, outdir, logObject)
	logObject.info("Done prepping FASTA inputs!")

	# Step 1: BLASTn query scaffolds to reference scaffold.
	logObject.info("Running BLASTn analysis.")
	blast_results_file = ConSequences.run_blastn(ref_fasta_file, query_multifasta_file, outdir, logObject)
	logObject.info("Successfully ran BLASTn analysis.")

	# Step 2: Perform sliding window analysis
	logObject.info("Running Sliding Window Analysis (SWA) to Identify Conserved and Contiguous Windows Along Reference.")
	swa_results_file = ConSequences.run_swa(blast_results_file, ref_fasta_file, scan_window_size, scan_slide_step, matches_percentage, outdir, logObject)
	logObject.info("Successfully ran SWA analysis.")

	# Step 3: Cluster windows and order along reference
	logObject.info("Starting clustering of windows into adjacent groupings.")
	window_clusters, loc_to_samples = ConSequences.cluster_windows(swa_results_file, scan_window_size, scan_slide_step, matches_percentage, logObject)
	logObject.info("Successfully clustered windows into adjacent groupings.")

	# Step 4: Delineate segments
	logObject.info("Delineating segments on sliding window anlaysis.")
	ConSequences.delineate_segments(window_clusters, loc_to_samples, ref_fasta_file, unfilter_segs, scan_slide_step, outdir, logObject)
	logObject.info("Successfully delineated segments from window clusters.")

	logObject.info("Program successfully ran! Exiting!")
	ConSequences.closeLoggerObject(logObject)
	sys.exit(0)

main()