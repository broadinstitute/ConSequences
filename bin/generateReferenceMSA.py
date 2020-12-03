#!/usr/bin/env python

### Program: generateReferenceMSA.py
### Author: Rauf Salamzade
### The Broad Institute of MIT and Harvard
### Earl Lab / Bacterial Genomics Group

"""
Special thanks to Terrance Shea and Bruce Walker for extensive feedback and suggestions.
"""

import os
import sys
import argparse
from ConSequences import ConSequences

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: generateReferenceMSA.py
	Author: Rauf Salamzade
	The Broad Institute of MIT and Harvard
	Earl Lab / Bacterial Genomics Group

	This program will generate a reference-based MSA for a segment of interest using results from delineateSegmentsOnReference.py. 
	If facing difficulties, please raise issues on the github page: https://github.com/broadinstitute/consequences
	""", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('-r', '--ref_fasta', type=str, help="FASTA for reference scaffold upon which segment lies.", required=True)
	parser.add_argument('-s', '--start_coord', type=int, help="Starting coordinate of segment.", required=True)
	parser.add_argument('-e', '--end_coord', type=int, help="Ending coordinate of segment.", required=True)
	parser.add_argument('-m', '--mapping_scaffs', type=str, help="List of scaffolds with segment. One per line. Or comma separated list provided in between quotes.", required=True)
	parser.add_argument('-w', '--sliding_window_results', type=str, help="Sliding window results file which contains variant information.", required=True)
	parser.add_argument('-o', '--msa_output', type=str, help="Multiple-sequence-alignment to be used for rapid identification of signature sequences.", required=True)
	parser.add_argument('-l', '--log_file', type=str, help="Path to logging output file", required=False, default="generate_msa.log")
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
	start_coord = myargs.start_coord
	end_coord = myargs.end_coord
	mapping_scaffs_input = myargs.mapping_scaffs
	swa_results_file = os.path.abspath(myargs.sliding_window_results)
	msa_output_file = os.path.abspath(myargs.msa_output)
	log_file = os.path.abspath(myargs.log_file)

	if os.path.isfile(msa_output_file):
		sys.stderr.write("Output multiple-sequence-alignment file exists. Overwriting in 5 seconds ...\n ")
		#sleep(5))

	try:
		assert(ConSequences.is_fasta(ref_fasta_file) and ConSequences.is_fasta(swa_results_file))
	except:
		raise RuntimeError("Reference FASTA or Sliding-Window Results file(s) are invalid. Raising error.")

	"""
	START WORKFLOW
	"""

	# create logging object
	logObject = ConSequences.createLoggerObject(log_file)

	# Step 0: Log input arguments and update reference and query FASTA files.

	logObject.info("Saving parameters for future provedance.")
	parameter_values = [ref_fasta_file, start_coord, end_coord, swa_results_file, msa_output_file]
	parameter_names = ["Reference FASTA file", "Start Coordinate of Segment", "End Coordinate of Segment", "Sliding Window Results file", "Multiple Sequence Alignment Output file"]
	ConSequences.logParametersToObject(logObject, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# Step 1: Read/load in data.
	logObject.info("Gathering segment sequence along reference and load scaffolds mapping to segment.")
	ref_scaff_len, ref_seg_seq = ConSequences.gather_sequence(ref_fasta_file, start_coord, end_coord, logObject)
	mapping_scaffs = [ms.strip() for ms in mapping_scaffs_input.split(',')]
	if os.path.isfile(mapping_scaffs_input):
		mapping_scaffs = ConSequences.load_mapping_scaffs(mapping_scaffs_input, logObject)
	mismatch_data = ConSequences.load_mismatch_data(swa_results_file, mapping_scaffs, logObject)
	logObject.info("Successfully loaded in data.")

	# Step 2: Create msa
	logObject.info("Create Multiple Sequence Alignment.")
	ConSequences.create_msa(ref_seg_seq, ref_scaff_len, start_coord, mapping_scaffs, mismatch_data, msa_output_file, logObject)
	logObject.info("Successfully ran SWA analysis.")

	logObject.info("Program successfully ran! Exiting!")
	ConSequences.closeLoggerObject(logObject)
	sys.exit(0)

main()