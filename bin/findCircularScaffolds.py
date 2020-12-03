import os
import sys
import argparse
import pysam
from Bio import SeqIO
import subprocess
from time import sleep
from collections import defaultdict
from ConSequences import ConSequences

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program which uses bowtie2 reflexive mapping of reads and blastn reflexive alignment of scaffolds to assess whether 
	any scaffolds in a sample's assembly are circular.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-s', '--sample_id', help='name of sample. "-" cannot be in name.', required=True)
	parser.add_argument('-a', '--assembly', help='Assembly (FASTA format).', required=True)
	parser.add_argument('-1', '--forward_reads', help='forward reads of paired-end set (FASTQ format).', required=True)
	parser.add_argument('-2', '--reverse_reads', help='reverse reads of paired-end set (FASTQ format).', required=True)
	parser.add_argument('-o', '--outdir', help='specify output directory', required=True)
	parser.add_argument('-t', '--threads', type=int, help='threads for bowtie2', required=False, default=1)
	args = parser.parse_args()
	return args

def main():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE ARGUMENTS
	"""

	myargs = create_parser()
	sample = myargs.sample_id
	assembly = myargs.assembly
	reads1 = myargs.forward_reads
	reads2 = myargs.reverse_reads
	threads = myargs.threads
	outdir = os.path.abspath(myargs.outdir) + '/'

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory exists. Overwriting in 5 seconds ...\n ")
		sleep(5)
	else:
		os.system('mkdir %s' % outdir)

	try:
		assert(ConSequences.is_fasta(assembly) and os.path.isfile(reads1) and os.path.isfile(reads2))
		assembly = os.path.abspath(assembly)
		reads1 = os.path.abspath(reads1)
		reads2 = os.path.abspath(reads2)
	except:
		raise RuntimeError("Input FASTA files are invalid. Raising error.")

	try:
		assert(not '-' in sample)
	except:
		raise RuntimeError('Sample had "-" in name, please change and retry!')

	"""
	START WORKFLOW
	"""

	# create logging object
	log_file = outdir + 'Progress.log'
	logObject = ConSequences.createLoggerObject(log_file)

	# Step 0: Log input arguments and update reference and query FASTA files.

	logObject.info("Saving parameters for future provedance.")
	parameters_file = outdir + 'Parameter_Inputs.txt'
	parameter_values = [assembly, reads1, reads2, outdir, threads]
	parameter_names = ["Assembly", "Forward Reads", "Reverse Reads", "Output Directory", "Threads"]
	ConSequences.logParametersToFile(parameters_file, parameter_names, parameter_values)
	logObject.info("Done saving parameters!")

	# Step 1: Perform Reflexive Alignment of Assembly with BLASTn
	logObject.info("Running BLASTn analysis to see if scaffold ends align and perform trimming.")
	circular_data, scaff_lens, trimmed_local_assembly = ConSequences.create_local_copy_and_trim_assembly(sample, assembly, outdir, logObject)
	logObject.info("Successfully ran BLASTn analysis.")

	# Step 2: Perform Reflexive Alignment of Paired-end Reads to Assembly
	logObject.info("Running Bowtie2 reflexive alignment analysis to determine presence of bridging paired-end reads.")
	scaff_bridging_reads = ConSequences.get_bridging(sample, scaff_lens, trimmed_local_assembly, reads1, reads2, outdir, threads, logObject)
	logObject.info("Successfully ran reflexive bowtie2 analysis.")

	"""
	Report Results!!!
	"""

	outf = open(outdir + 'Circular_Scaffolds_Inference.txt', 'w')
	outf.write('sample\tscaffold\tscaffold_length_before_trimming\ttrimming_length\tmulti_side_mapping_read_support\n')
	for s in scaff_lens:
		printlist = [sample, s, scaff_lens[s], circular_data[s], scaff_bridging_reads[s]]
		printlist = [str(x) for x in printlist]
		outf.write('\t'.join(printlist) + '\n')
	outf.close()

main()