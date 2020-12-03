import os
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter
import logging
import subprocess
import math
import itertools
import binascii
import pysam

def create_local_copy_and_trim_assembly(sample, assembly, outdir, logObject):
	"""
	Function for findCircularScaffolds.py to create a local instance of the assembly and perform
	BLASTn analysis to see if scaffold ends align, perform trimming where needed if so.
	"""

	logObject.info("Parsing assembly:")
	seqlens = {}
	with open(assembly) as oa:
		for rec in SeqIO.parse(oa, 'fasta'):
			seqlens[rec.id] = len(str(rec.seq))
			logObject.info("Length of scaffold %s is %d" % (rec.id, len(str(rec.seq))))

	""" Run BLASTn Analysis + Trimming """

	try:
		blastn_result = outdir + 'self_blast.txt'
		blastn_cmd = 'blastn -subject %s -query %s -outfmt \"6 std qcovhsp\" -ungapped -max_target_seqs 100000 -out %s' % (assembly, assembly, blastn_result)
		logObject.info("Running following BLASTn command: %s" % blastn_cmd)
		subprocess.call(blastn_cmd, shell=True)
	except:
		logObject.error("Error running BLASTn! Exiting now ...")
		raise RuntimeError("Error running BLASTn! Exiting now ...")

	### Get self eating snakes
	circular_data = defaultdict(int)
	with open(blastn_result):
		for line in open(blastn_result):
			line = line.strip()
			ls = line.split('\t')
			if ls[0] == ls[1] and not (ls[-4] == ls[-6] and ls[-5] == ls[-7]) and not int(ls[-1]) == 100:
				if float(ls[2]) >= 99.0 and float(ls[-3]) <= 1e-5:
					qs, qe = [int(ls[-7]), int(ls[-6])]
					ss, se = [int(ls[-5]), int(ls[-4])]
					if (ss == 1 and qe == seqlens[ls[0]]) or (se == 1 and qs == seqlens[ls[0]]):
						cut = max(ss, se)
						min_que = min(qs, qe)
						if min_que < cut:
							cut = min_que
						if cut > circular_data[ls[0]]:
							circular_data[sample + '|' + ls[0]] = cut

	trimmed_assembly = outdir + sample + '.trimmed.fasta'
	trimmed_assembly_handle = open(trimmed_assembly, 'w')
	seqlens_trim = {}
	with open(assembly) as oa:
		for rec in SeqIO.parse(assembly, 'fasta'):
			if rec.id in circular_data:
				logObject.info("Trimming %d bp from the beginning of scaffold %s" % (circular_data[sample + '|' + rec.id], rec.id))
			trim_seq = str(rec.seq)[circular_data[sample + '|' + rec.id]:]
			trimmed_assembly_handle.write('>' + sample + '|' + rec.id + '\n' + trim_seq + '\n')
			seqlens_trim[sample + '|' + rec.id] = len(trim_seq)
	trimmed_assembly_handle.close()
	return([circular_data, seqlens_trim, trimmed_assembly])

def get_bridging(sample, scaff_lens, assembly, reads1, reads2, outdir, threads, logObject):
	"""
	Function to perform reflexive alignment of paired-end fragment reads to sample assembly and determine if there are
	bridging read pairs for any of the scaffolds.
	"""

	""" Run Bowtie2 Analysis """
	bowtie2_index = outdir + 'reference'
	sam_file = outdir + 'mapped.sam'
	bam_file = outdir + 'mapped.bam'
	bam_file_sorted = outdir + 'mapped.sorted.bam'
	multisided_bam_file_sorted = outdir + 'bridging.sorted.bam'

	try:
		bowtie2_build_cmd = ['bowtie2-build', assembly, bowtie2_index]
		logObject.info("Running following bowtie2 command: %s" % " ".join(bowtie2_build_cmd))
		subprocess.call(' '.join(bowtie2_build_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

		bowtie2_cmd = ['bowtie2', '--very-sensitive-local', '-a', '-x', bowtie2_index, '-1', reads1, '-2', reads2, '-S',
					   sam_file, '-p', str(threads)]
		logObject.info("Running following bowtie2 command: %s" % " ".join(bowtie2_cmd))
		subprocess.call(' '.join(bowtie2_cmd), shell=True)

		sam_to_bam_cmd = ["samtools", "view", "-Sb", sam_file]
		ob = open(bam_file, 'w')
		logObject.info("Running following samtools command: %s" % " ".join(sam_to_bam_cmd))
		subprocess.call(' '.join(sam_to_bam_cmd), stdout=ob, shell=True)
		ob.close()

		sort_bam_cmd = ["samtools", "sort", '-@', str(threads), bam_file, '-o', bam_file_sorted]
		logObject.info("Running following samtools command: %s" % " ".join(sort_bam_cmd))
		subprocess.call(' '.join(sort_bam_cmd), shell=True)

		index_bam_cmd = ["samtools", "index", bam_file_sorted]
		logObject.info("Running following samtools command: %s" % " ".join(index_bam_cmd))
		subprocess.call(' '.join(index_bam_cmd), shell=True)
	except:
		logObject.error("Error running Bowtie2 or downstream processing with samtools! Exiting now ...")
		raise RuntimeError("Error running Bowtie2 or downstream processing with samtools! Exiting now ...")


	scaff_multi_side_reads = defaultdict(int)
	try:
		bam_handle = pysam.AlignmentFile(bam_file_sorted, 'rb')
		multisided_bam = outdir + "multisided.bam"
		pairedreads = pysam.AlignmentFile(multisided_bam, "wb", template=bam_handle)
		for s, slen in scaff_lens.items():
			for read1s, read2s in multi_sided_read_pair_identifier(bam_handle, region_string=s, slen=slen):
				scaff_multi_side_reads[s] += 1
				for r1 in read1s:
					pairedreads.write(r1)
				for r2 in read2s:
					pairedreads.write(r2)
		pairedreads.close()
		bam_handle.close()
	except:
		logObject.info("Error parsing briding paired-end reads! Exiting now ...")
		raise RuntimeError("Error parsing briding paired-end reads! Exiting now ...")

	try:
		sort_bam_cmd = ["samtools", "sort", '-@', str(threads), multisided_bam, '-o', multisided_bam_file_sorted]
		logObject.info("Running following samtools command: %s" % " ".join(sort_bam_cmd))
		subprocess.call(' '.join(sort_bam_cmd), shell=True)

		index_bam_cmd = ["samtools", "index", multisided_bam_file_sorted]
		logObject.info("Running following samtools command: %s" % " ".join(index_bam_cmd))
		subprocess.call(' '.join(index_bam_cmd), shell=True)
	except:
		logObject.error("Error running sorting, indexing of filtered bam file containing just the bridging paired-end reads! Exiting now ...")
		raise RuntimeError("Error running sorting, indexing of filtered bam file containing just the bridging paired-end reads! Exiting now ...")

	try:
		cleanup_cmd = ['rm', '-f', sam_file, bam_file, bam_file_sorted]
		logObject.info("Running cleanup of intermediate bam files!")
		subprocess.call(cleanup_cmd)
	except:
		logObject.error("Error running cleanup! Exiting now ...")
		raise RuntimeError("Error running cleanup! Exiting now ...")

	return scaff_multi_side_reads

def multi_sided_read_pair_identifier(bam, region_string=None, slen=None):
	"""
	Function which reads through bam and identifies if a read pair supports multi-sided mapping/bridging of scaffold ends.
	"""
	read_map_counts = defaultdict(int)

	read_dict = defaultdict(lambda: defaultdict(list))
	for i, read in enumerate(bam.fetch(region=region_string)):
		readlen = read.reference_length
		if not readlen: continue
		qname = read.query_name

		location = read.reference_start + 1
		location = int(location)
		if read.is_read1:
			read_map_counts[qname] += 1
		insertion_location_left = location
		insertion_location_right = location + readlen
		if abs(insertion_location_right - slen) < 500:
			if read.is_read1:
				read_dict[qname][0].append(['end', read])
			else:
				read_dict[qname][1].append(['end', read])
		elif abs(insertion_location_left - 1) < 500:
			if read.is_read1:
				read_dict[qname][0].append(['start', read])
			else:
				read_dict[qname][1].append(['start', read])

	results = []
	for qname in read_dict:
		r1_reads = [x[1] for x in read_dict[qname][0]]
		r2_reads = [x[1] for x in read_dict[qname][1]]
		r1_sides = set([x[0] for x in read_dict[qname][0]])
		r2_sides = set([x[0] for x in read_dict[qname][1]])
		if ('start' in r1_sides and 'end' in r2_sides) or ('end' in r1_sides and 'start' in r2_sides):
			results.append(tuple([r1_reads, r2_reads]))
	return results

def run_prep(ref_fasta_file, query_fasta_file, scan_window_size, ref_circular, query_circular, outdir, logObject):
	"""
	Update input fasta files to account for circularity where possible.

	:param ref_fasta_file: reference fasta input (un-edited)
	:param query_fasta_file: query fasta input (un-edited)
	:param ref_circular: boolean for whether reference is circular
	:param query_circular: set of query scaffold identifiers representing set which are circular
	:param outdir: primary output directory path
	:param logObject: logging logger object
	:return: [updated reference fasta file path, updated query fasta file path]
	"""

	step_dir = outdir + "Processed_Input_FASTA/"
	if not os.path.isdir(step_dir): os.system('mkdir %s' % step_dir)

	proc_ref_fasta_file = step_dir + 'reference_scaffold.fasta'
	proc_query_fasta_file = step_dir + 'query_scaffolds.fasta'

	proc_ref_fasta_handle = open(proc_ref_fasta_file, 'w')
	num_seqs = 0
	with open(ref_fasta_file) as orff:
		for rec in SeqIO.parse(orff, 'fasta'):
			seq = str(rec.seq)
			seqlen = len(seq)
			try:
				assert(seqlen >= scan_window_size)
			except:
				logObject.error("Reference scaffold is shorter than window size to be used in sliding window analysis. Exiting now...")
				raise RuntimeError("Reference scaffold is shorter than window size to be used in sliding window analysis. Exiting now...")
			proc_ref_fasta_handle.write('>' + rec.id + '|' + str(seqlen) + '\n')
			if ref_circular:
				logObject.info("Reference scaffold %s is regarded as complete and circular." % rec.id)
				proc_ref_fasta_handle.write(seq + seq[:scan_window_size] + '\n')
			else:
				logObject.info("Reference scaffold %s is not regarded as complete and circular." % rec.id)
				proc_ref_fasta_handle.write(seq + '\n')
			num_seqs += 1
	proc_ref_fasta_handle.close()

	try:
		assert(num_seqs == 1)
	except:
		logObject.error("More than one scaffold sequence in reference FASTA file not allowed. Exiting now...")
		raise RuntimeError("More than one scaffold sequence in reference FASTA file not allowed. Exiting now...")

	proc_query_fasta_handle = open(proc_query_fasta_file, 'w')
	with open(query_fasta_file) as oqff:
		for rec in SeqIO.parse(oqff, 'fasta'):
			seq = str(rec.seq)
			seqlen = len(seq)
			proc_query_fasta_handle.write('>' + rec.id + '|' + str(seqlen) + '\n')
			if seqlen < scan_window_size:
				logObject.warning("Query scaffold %s is shorter than window size to be used in sliding window analysis and dropped from analysis." % (rec.id))
				continue
			if rec.id in query_circular:
				logObject.info("Query scaffold %s is regarded as complete and circular." % rec.id)
				proc_query_fasta_handle.write(seq + seq[:scan_window_size] + '\n')
			else:
				logObject.info("Query scaffold %s is not regarded as complete and circular." % rec.id)
				proc_query_fasta_handle.write(seq + '\n')
	proc_query_fasta_handle.close()

	try:
		assert(is_fasta(proc_ref_fasta_file) and is_fasta(proc_query_fasta_file))
	except:
		logObject.error("Processing of input FASTAs had issues. Exiting now...")
		raise RuntimeError("Processing of input FASTAs had issues. Exiting now...")

	return([proc_ref_fasta_file, proc_query_fasta_file])

def run_blastn(reference_fasta, query_fasta, outdir, logObject):
	"""
	Run BLASTn analysis.

	:param reference_fasta: reference fasta path
	:param query_fasta: query fasta path
	:param outdir: primary output directory path
	:param logObject: logging logger object
	:return: path to BLASTn output in output format 5.
	"""

	blastn_results = outdir + 'BLAST_results.txt'

	blastn_cmd = ['blastn', '-subject', reference_fasta, '-query', query_fasta, '-evalue', '1e-5',
				  '-max_target_seqs', '100000', '-outfmt', '5', '-out', blastn_results]


	logObject.info("Executing the command: %s" % ' '.join(blastn_cmd))
	stdout_file = outdir + 'blast.log';	stderr_file = outdir + 'blast.err'
	stdout_handle = open(stdout_file, 'w'); stderr_handle = open(stderr_file, 'w')
	proc = subprocess.Popen(blastn_cmd, stdout=stdout_handle, stderr=stderr_handle)
	proc.wait()
	stdout_handle.close(); stderr_handle.close()
	try:
		assert(os.path.isfile(blastn_results) and os.path.getsize(blastn_results) > 100)
	except:
		logObject.error("Error running BLASTn analysis, check logs and input FASTA. Exiting now...")
		raise RuntimeError("Error running BLASTn analysis, check logs and input FASTA. Exiting now...")

	return(blastn_results)

def run_swa(blast_results_file, reference_fasta, scan_window_size, scan_slide_step, matches_percentage, outdir, logObject):
	"""
	Run sliding window analysis to identify conserved windows.
	:param blast_results_file: BLASTn xml results file.
	:param reference_fasta: The reference fasta file.
	:param scan_window_size: The sliding window size.
	:param scan_slide_step: The step size for sliding window analysis.
	:param matches_percentage: The identity required for reference window to be considered present in query scaffold.
	:param outdir: path primary output directory
	:param logObject: logging logger object.
	:return: path to SWA results file.
	"""
	logObject.info("Parsing BLASTn XML results file.")
	matches_cutoff = int(math.floor((matches_percentage/100.0)*scan_window_size))
	sub_hits = []

	""" Parse BLASTn XML results file to extract valid HSPs for consideration."""
	with open(blast_results_file) as obxf:
		query_name = None
		hsp_data = {}
		for line in obxf:
			line = line.strip()
			if line.startswith("<Iteration_query-def>"):
				query_name = line.split('>')[1].split('<')[0]
			elif line.startswith("<Hsp_query-from>"):
				hsp_data['qstart'] = int(line.split('>')[1].split('<')[0])
			elif line.startswith("<Hsp_query-to>"):
				hsp_data['qend'] = int(line.split('>')[1].split('<')[0])
			elif line.startswith("<Hsp_hit-from>"):
				hsp_data['sstart'] = int(line.split('>')[1].split('<')[0])
			elif line.startswith("<Hsp_hit-to>"):
				hsp_data['send'] = int(line.split('>')[1].split('<')[0])
			elif line.startswith("<Hsp_identity>"):
				hsp_data['identity'] = int(line.split('>')[1].split('<')[0])
			elif line.startswith("<Hsp_qseq>"):
				hsp_data['query'] = line.split('>')[1].split('<')[0]
			elif line.startswith("<Hsp_hseq>"):
				hsp_data['subject'] = line.split('>')[1].split('<')[0]
			elif line.startswith("<Hsp_midline>"):
				hsp_data['matches'] = line.split('>')[1].split('<')[0]
			elif line.startswith("</Hsp>"):
				if hsp_data['identity'] >= matches_cutoff:
					subseq = hsp_data['subject']
					queseq = hsp_data['query']
					matseq = hsp_data['matches']
					match = ""
					query = ""
					subject = ""
					insertions = defaultdict(int)
					insertion_snip = ""
					ref_pos = 0
					for p, bp in enumerate(subseq):
						if bp != '-':
							if len(insertion_snip) != 0:
								if hsp_data['sstart'] < hsp_data['send']:
									insertions[tuple([ref_pos + hsp_data['sstart'] - 1,
													  ref_pos + hsp_data['sstart']])] = insertion_snip
								else:
									insertions[tuple(
										[hsp_data['sstart'] - ref_pos + 1, hsp_data['sstart'] - ref_pos])] = revcomp(insertion_snip)
								insertion_snip = ""
							match += matseq[p]
							query += queseq[p]
							subject += bp
							ref_pos += 1
						else:
							insertion_snip += queseq[p]

					sub_hits.append(
						[query_name, hsp_data['sstart'], hsp_data['send'], hsp_data['qstart'], hsp_data['qend'], match,
						 query, subject, insertions])
				hsp_data = {}

	logObject.info("Processed - %d HSPs in consideration." % len(sub_hits))

	logObject.info("Starting sliding window analysis (SWA)!")
	total_windows = 0
	self_mapping_windows = 0
	multi_mapping_windows = 0
	swa_results_file = outdir + 'SWA_Results.txt'
	swa_results_handle = open(swa_results_file, 'w')
	swa_results_handle.write('\t'.join(
		['# reference scaffold ID', 'window start', 'window stop', 'number of matching query scaffolds',
		'matching query scaffolds', 'mismatch information for matching matching query scaffolds']))

	""" Run Sliding Window Analysis """
	try:
		with open(reference_fasta) as osff:
			for rec in SeqIO.parse(osff, 'fasta'):
				reference_id = rec.id
				seqlen = len(str(rec.seq))
				start = 1
				stop = scan_window_size
				while stop <= seqlen:
					total_windows += 1
					matches = set([])
					mismatches_for_matches = set([])

					total_overlap = defaultdict(int)
					total_overlap_coords = defaultdict(set)
					isl_range = set(range(start, stop + 1))

					for hsp in sorted(sub_hits, key=itemgetter(6), reverse=True):
						que, hit_start, hit_stop, que_start, que_stop, pos_matches, pos_query, pos_subject, pos_insertions = hsp
						if stop < min(hit_start, hit_stop) or max(hit_start, hit_stop) < start: continue
						max_possible_identities = len(
							set(range(min(hit_start, hit_stop), max(hit_start, hit_stop))).intersection(isl_range))
						if max_possible_identities < matches_cutoff: continue

						direction = '+'
						if hit_start > hit_stop: pos_matches = pos_matches[::-1]; pos_subject = revcomp(pos_subject); direction = '-'; pos_query = revcomp(pos_query);

						hit_start, hit_stop = sorted([hit_start, hit_stop])

						insertion_score = 0
						insertion_coords = set([])
						for ips in pos_insertions.items():
							if ips[0][0] in isl_range and ips[0][1] in isl_range:
								insertion_score += len(ips[1])
								insertion_coords.add('%d.%d.%s' % (ips[0][0], ips[0][1], ips[1]))

						overlap = 0
						hit_range = set([])
						mismatch = set([])
						mis_range = set([])

						for p, pos in enumerate(range(hit_start, hit_stop + 1)):
							# print(str(pos) + '\t' + pos_matches[p] + '\t' + pos_subject[p] + '\t' + pos_query[p])
							if pos_matches[p] == '|':
								hit_range.add(pos)
							elif pos in isl_range:
								mis_range.add(pos)
								mismatch.add(str(pos) + '_' + pos_query[p] + '_' + pos_subject[p] + '_' + direction)
						if insertion_score > 0: mismatch.add("INSERTION_" + '_'.join(sorted(insertion_coords)))

						overlap = max(len(isl_range.intersection(hit_range)) - insertion_score, 0)

						missing_positions = isl_range.difference(hit_range).difference(mis_range)
						for pos in missing_positions:
							mismatch.add(str(pos) + '_-_?_?')

						if total_overlap[que] < overlap:
							total_overlap[que] = overlap
							total_overlap_coords[que] = mismatch

					for q in total_overlap.keys():
						island_coverage_total = total_overlap[q]
						if island_coverage_total >= matches_cutoff:
							matches.add(q)
							mismatches_for_matches.add(
								q + ': ' + ','.join([str(x) for x in sorted(total_overlap_coords[q])]))

					if len(matches) >= 1:
						try:
							assert (reference_id in matches)
							self_mapping_windows += 1
						except:
							start += scan_slide_step
							stop += scan_slide_step
							continue
					if len(matches) >= 2:
						try:
							assert (reference_id in matches)
						except:
							logObject.error("No self mapping hit yet another scaffold maps to this window. Not supposed to happen. Exiting now...")
							raise RuntimeError("No self mapping hit yet another scaffold maps to this window. Not supposed to happen. Exiting now...")

						multi_mapping_windows += 1
						swa_results_handle.write('\t'.join([str(x) for x in [reference_id, start, stop, len(matches), '; '.join(sorted(matches)), '; '.join(sorted(mismatches_for_matches))]]) + '\n')
					start += scan_slide_step
					stop += scan_slide_step
	except:
		logObject.error("Had an issue when running sliding window analysis. Exiting now...")
		raise RuntimeError("Had an issue when running sliding window analysis. Exiting now...")

	swa_results_handle.close()
	logObject.info('SWA Results - Total windows: %d' % total_windows)
	logObject.info('SWA Results - Self-mapping windows: %d' % self_mapping_windows)
	logObject.info('SWA Results - Multi-mapping windows: %d' % multi_mapping_windows)

	try:
		assert(multi_mapping_windows > 0)
	except:
		logObject.info("No conserved windows identified upon reference scaffold with query scaffolds provided.")
		logObject.info("Program successfully ran! Exiting!")
		sys.exit(0)

	return swa_results_file

def cluster_windows(swa_results_file, scan_window_size, scan_slide_step, matches_percentage, logObject):
	"""
	Cluster/stitch adjacent windows into groupings.
	:param swa_results_file: Results file from SWA.
	:param outdir: path to primary output directory.
	:param logObject: logging logger object.
	:return: [ list of lists containing single-linkage clusters of adjacent windows, loc_to_samples]
	"""
	matches_cutoff = int(math.floor((matches_percentage/100.0)*scan_window_size))

	try:
		""" Read in SWA results """
		reflen = None
		ref_windows = []
		loc_to_samples = {}
		overlap_pairs = []
		with open(swa_results_file) as osrf:
			for line in osrf:
				if line.startswith("#"): continue
				line = line.strip()
				ls = line.split('\t')
				reflen = int(ls[0].split('|')[-1])

				start = int(ls[1])
				stop = int(ls[2])

				if start >= reflen: continue
				if stop > reflen:
					stop = stop - reflen

				loc_key = '%d|%d' % (start, stop)
				locs = [start, stop]
				samples = set(ls[-2].split('; '))
				ref_windows.append([locs, samples])
				loc_to_samples[loc_key] = samples
				overlap_pairs.append([loc_key, loc_key])

		"""	Cluster overlapping windows	"""

		for i, w1 in enumerate(ref_windows):
			for j, w2 in enumerate(ref_windows):
				if i < j:
					w1_range = None;
					w2_range = None
					if w1[0][0] >= w1[0][1]: w1_range = set(range(w1[0][0], reflen + 1)).union(set(range(1, w1[0][1] + 1)))
					if w2[0][0] >= w2[0][1]: w2_range = set(range(w2[0][0], reflen + 1)).union(set(range(1, w2[0][1] + 1)))

					w1key = '%d|%d' % (w1[0][0], w1[0][1])
					w2key = '%d|%d' % (w2[0][0], w2[0][1])

					if w1_range or w2_range:
						if w1_range and not w2_range:
							w2_range = set(range(w2[0][0], w2[0][1] + 1))
						if w2_range and not w1_range:
							w1_range = set(range(w1[0][0], w1[0][1] + 1))
						if len(w1_range.intersection(w2_range)) >= matches_cutoff:
							overlap_pairs.append(sorted([w1key, w2key]))
					elif (w1[0][0] + scan_slide_step) == w2[0][0]:
							overlap_pairs.append(sorted([w1key, w2key]))
		"""	
		Solution for single-linkage clustering taken from mimomu's repsonse in the stackoverflow page:
		https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements?lq=1
		"""
		L = overlap_pairs
		LL = set(itertools.chain.from_iterable(L))
		for each in LL:
			components = [x for x in L if each in x]
			for i in components:
				L.remove(i)
			L += [list(set(itertools.chain.from_iterable(components)))]
		return [L, loc_to_samples]

	except:
		logObject.error("Had problem when clustering windows. This is likely an issue with the program. Please raise github issue and share example input and results from SWA reslts file. Now exiting ...")
		raise RuntimeError("Had problem when clustering windows. This is likely an issue with the program. Please raise github issue and share example input and results from SWA reslts file. Now exiting ...")

def delineate_segments(L, loc_to_samples, reference_fasta, unfilter_segs, scan_slide_step, outdir, logObject):
	"""
	Delineate segments within each cluster of windows.
	:param L: single-linkage clustering of windows.
	:param loc_to_samples: window id to samples set dictionary.
	:param reference_fasta: reference fasta
	:param unfilter_segs: whether to filter or segments to make more concise or not.
	:param outdir: path to primary output directory.
	:param logObject: logging logger object.
	:return:
	"""
	reflen = None
	reference_id = None
	with open(reference_fasta) as osff:
		for rec in SeqIO.parse(osff, 'fasta'):
			reference_id = rec.id
			reflen = len(str(rec.seq))

	segsum_file = open(outdir + 'Segment_Results.txt', 'w')
	elements_file = open(outdir + 'Segment_Windows.txt', 'w')
	clusters_file = open(outdir + 'Clustered_Windows.txt', 'w')
	for ci, c in enumerate(L):
		bridging_el = [[int(x.split('|')[0]), int(x.split('|')[1])] for x in c if
					   int(x.split('|')[0]) > int(x.split('|')[1])]
		non_bridging_el = [[int(x.split('|')[0]), int(x.split('|')[1])] for x in c if
						   int(x.split('|')[0]) < int(x.split('|')[1])]

		# create flat list of endpoints of windows in cluster
		nb_vals_flat = list(itertools.chain(*non_bridging_el)) + list(itertools.chain(*bridging_el))
		nb_vals_flat = [x for x in nb_vals_flat if x <= reflen] + [(x - reflen) for x in nb_vals_flat if x > reflen]

		left_el = []; right_el = []
		if len(bridging_el) > 0:

			cutoffs = []
			cutoff_met = 0
			prev_val = None
			for val in sorted(nb_vals_flat):
				if prev_val:
					if abs(val - prev_val) >= scan_slide_step+1:
						cutoff_met += 1
						cutoffs.append(val - 1)
				prev_val = val

			try:
				assert (cutoff_met <= 1)
			except:
				try:
					assert (cutoff_met <= 2)
				except:
					logObject.error("Should not happen. Raise github issue and share input/SWA results.")
					raise RuntimeError("Should not happen. Raise github issue and share input/SWA results.")

			if len(cutoffs) > 0:
				cutoff_min = sorted(cutoffs)[0]
				left_el = [x for x in non_bridging_el if x[0] > cutoff_min]
				right_el = [x for x in non_bridging_el if x[0] < cutoff_min]
			else:
				left_el = non_bridging_el
		else:
			left_el = non_bridging_el

		window_index = 1
		windex_to_loc = {}

		for w in sorted(left_el, key=lambda x: x[0]):
			windex_to_loc[window_index] = '%d|%d' % (w[0], w[1])
			window_index += 1
		for w in sorted(bridging_el, key=lambda x: x[0]):
			windex_to_loc[window_index] = '%d|%d' % (w[0], w[1])
			window_index += 1
		for w in sorted(right_el, key=lambda x: x[0]):
			windex_to_loc[window_index] = '%d|%d' % (w[0], w[1])
			window_index += 1

		all_cluster_windows = []
		segment_windows = defaultdict(list)
		for i, w1 in enumerate(sorted(windex_to_loc.keys())):
			wind1 = windex_to_loc[w1]
			all_cluster_windows.append(reference_id + '|' + wind1)
			w1_samples = set(loc_to_samples[wind1])
			curr_intersection_set = w1_samples
			for j, w2 in enumerate(sorted(windex_to_loc.keys())):
				if i <= j:
					wind2 = windex_to_loc[w2]
					w2_samples = set(loc_to_samples[wind2])
					curr_intersection_set = curr_intersection_set.intersection(w2_samples)
					if len(curr_intersection_set) > 0:
						segment_windows[tuple(sorted(list(curr_intersection_set)))].append(sorted([w1, w2]))
					else:
						break

		for cs in segment_windows:
			cs_len = len(set(cs))
			overlap_windows = segment_windows[cs]
			EL = overlap_windows
			ELL = set(itertools.chain.from_iterable(EL))
			for each in ELL:
				components = [x for x in EL if each in x]
				for i in components:
					EL.remove(i)
				EL += [list(set(itertools.chain.from_iterable(components)))]
			for si, seg in enumerate(EL):
				windows = []
				matches_cs =  False
				for wi in sorted(seg):
					windows.append(windex_to_loc[wi])
					wi_samples = set(loc_to_samples[windex_to_loc[wi]])
					if len(set(wi_samples).union(set(cs))) == cs_len:
						matches_cs = True
				if (matches_cs and not unfilter_segs) or unfilter_segs:
					elements_file.write(reference_id + '\t' + str(ci + 1) + '\t' + str(si + 1) + '\t' + ', '.join(windows) + '\t' + ', '.join(sorted(list(cs))) + '\n')

					start = None
					stop = None
					covered = set([])
					for i, w in enumerate(windows):
						ws = [int(x) for x in w.split('|')]
						if i == 0: start = ws[0]
						stop = ws[1]
						if ws[0] < ws[1]:
							for p in range(ws[0], ws[1] + 1):
								covered.add(p)
						else:
							for p in range(ws[0], reflen + 1):
								covered.add(p)
							for p in range(1, ws[1] + 1):
								covered.add(p)

					if reflen == len(covered):
						stop = reflen

					segsum_file.write('|'.join(reference_id.split('|')[:-1]) + '\t' + str(start) + '\t' + str(stop) + '\t' + ', '.join(sorted(list(cs))) + '\n')

		clusters_file.write(' '.join(all_cluster_windows) + '\n')

	elements_file.close()
	clusters_file.close()
	segsum_file.close()

def gather_sequence(ref_fasta, start, stop, logObject):
	seq = ""
	with open(ref_fasta) as orf:
		for rec in SeqIO.parse(orf, 'fasta'):
			seq = str(rec.seq)
	seqlen = len(seq)
	subseq = ''
	if max([start, stop]) > seqlen:
		logObject.error('Start/end coordinates don\'t seem to align with lenght of reference fasta. Exiting now.')
		raise RuntimeError('Start/end coordinates don\'t seem to align with lenght of reference fasta. Exiting now.')

	if start > stop:
		pos_range = set(range(start - 1, seqlen))
		for p, bp in enumerate(seq):
			if p in pos_range: subseq += bp
		pos_range = set(range(0, stop))
		for p, bp in enumerate(seq):
			if p in pos_range: subseq += bp
	else:
		pos_range = set(range(start - 1, stop))
		for p, bp in enumerate(seq):
			if p in pos_range: subseq += bp
	return([len(seq), subseq])

def load_mapping_scaffs(mapping_scaffs_file, logObject):
	mapping_scaffs = set([])
	with open(mapping_scaffs_file) as omsf:
		for line in omsf:
			line = line.strip()
			mapping_scaffs.add(line)
	return(mapping_scaffs)

def load_mismatch_data(swa_results_file, mapping_scaffs, logObject):
	logObject.info("Loading mismatch data from sliding window results file.")
	mismatch_data = defaultdict(lambda: defaultdict(dict))
	with open(swa_results_file) as owf:
		for l in owf:
			l = l.strip()
			if l.startswith('#'): continue
			ls = l.split('\t')
			mismatches = ls[-1].split('; ')
			for sm in mismatches:
				scaff = sm.split(':')[0]
				coding_info = sm.split(':')[1].strip().split(',')
				if scaff in mapping_scaffs and coding_info[0].strip() != '':
					for mm in coding_info:
						if not mm.startswith("INSERTION"):
							pos,alt,ref,dir = mm.split('_')
							pos = int(pos)
							mismatch_data[scaff]['reference'][pos] = ref
							mismatch_data[scaff]['alternate'][pos] = alt
							mismatch_data[scaff]['direction'][pos] = dir
						else:
							for ins in mm.split('_')[1:]:
								ins_pos1, ins_pos2, ins_len =  ins.split('.')
								ins_pos1 = int(ins_pos1)
								ins_pos2 = int(ins_pos2)
								ins_len = len(ins_len)
								mismatch_data[scaff]['reference'][ins_pos1] = '?'
								mismatch_data[scaff]['reference'][ins_pos2] = '?'
								mismatch_data[scaff]['direction'][ins_pos1] = '?'
								mismatch_data[scaff]['direction'][ins_pos2] = '?'
								mismatch_data[scaff]['alternate'][ins_pos1] = '-'
								mismatch_data[scaff]['alternate'][ins_pos2] = '-'
	return mismatch_data

def create_msa(ref_seg_seq, ref_scaff_len, start_coord, mapping_scaffs, mismatch_data, msa_output_file, logObject):
		logObject.info("Starting to create multiple sequence alignment!")
		try:
			data = [sorted(mapping_scaffs)]
			for p, ref_bp in enumerate(ref_seg_seq):
				position_values = []
				pos = p+start_coord
				if pos > ref_scaff_len:
					pos = (p+start_coord)-ref_scaff_len
				for j, s in enumerate(sorted(mapping_scaffs)):
					if pos in mismatch_data[s]['reference'].keys():
						assert(ref_bp == mismatch_data[s]['reference'][pos] or mismatch_data[s]['reference'][pos] == '?')
						position_values.append(mismatch_data[s]['alternate'][pos])
					else:
						position_values.append(ref_bp)
				data.append(position_values)
		except:
			logObject.error("Problem with creating msa!")
			raise RuntimeError("Problem with creating msa!")
		logObject.info("Writing multiple sequence alignment to file.")
		try:
			result_msa_file = open(msa_output_file, 'w')
			for ls in zip(*data):
				result_msa_file.write('>' + str(ls[0]) + '\n' + ''.join(ls[1:]) + '\n')
			result_msa_file.close()
		except:
			logObject.error("Problem with writing msa to file!")
			raise RuntimeError("Problem with writing msa to file!")

def createLoggerObject(log_file):
	"""
	Function which creates logger object.
	:param log_file: path to log file.
	:return: logging logger object.
	"""
	logger = logging.getLogger('task_logger')
	logger.setLevel(logging.DEBUG)
	# create file handler which logs even debug messages
	fh = logging.FileHandler(log_file)
	fh.setLevel(logging.DEBUG)
	# create console handler with a higher log level
	ch = logging.StreamHandler()
	ch.setLevel(logging.ERROR)
	# create formatter and add it to the handlers
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M")
	fh.setFormatter(formatter)
	ch.setFormatter(formatter)
	# add the handlers to the logger
	logger.addHandler(fh)
	logger.addHandler(ch)
	#logger.handlers[0].stream = sys.stderr
	return logger

def closeLoggerObject(logObject):
	"""
	Function which closes/terminates loggerObject.
	:param logObject: logging logger object to close
	"""
	handlers = logObject.handlers[:]
	for handler in handlers:
		handler.close()
		logObject.removeHandler(handler)

def comp(seqx):
	seqxrc = Seq(seqx).complement()
	return (str(seqxrc))

def check_if_sample_has_segment(ref_msa, ref_seq, motif_msa, sample_kmers, kmer_length, output_prefix):
	try:
		start = 0
		stop = kmer_length
		mid_step_missed_count = 0
		tot_step_missed_count = 0
		tot_steps = 0
		while stop <= len(str(ref_seq)):
			gap_flag = False
			options = set([])
			for s in motif_msa:
				ss = s[start:stop]
				if 'N' in ss or '-' in ss: gap_flag = True; break
				options.add(ss)

			if gap_flag:
				start += 1
				stop += 1
				continue

			step_found = False
			for op in sorted(options):
				if op in sample_kmers:
					step_found = True
			if not step_found:
				if start > 100 and len(ref_seq) - stop > 100:
					mid_step_missed_count += 1
				tot_step_missed_count += 1
			tot_steps += 1
			start += 1
			stop += 1

		output_handle = open(output_prefix + '_results.txt', 'a+')
		output_handle.write("%s\t%d\t%d\t%d\n" % (ref_msa, tot_step_missed_count, mid_step_missed_count, tot_steps))
		output_handle.close()
	except:
		raise RuntimeError("ERROR: Problem with looking for segment assembly path in sample's reads/kmers.")

def kmerize_msa(ref_msa_file, reference, kmer_length):
	try:
		ref_seq = None
		motif_msa = []
		kmers = set([])
		with open(ref_msa_file) as omm:
			for rec in SeqIO.parse(omm, 'fasta'):
				start = 0
				stop = kmer_length
				while stop <= len(str(rec.seq)):
					kmer = str(rec.seq)[start:stop]
					kmers.add(kmer)
					kmers.add(revcomp(kmer))
					start += 1
					stop += 1
				if rec.id == reference:
					ref_seq = str(rec.seq)
				motif_msa.append(str(rec.seq))
		return ([ref_seq, motif_msa, kmers])
	except:
		raise RuntimeError('ERROR: Problem gather k-mer space for MSA!')

def kmerize_reads(illumina_read_files, output_prefix, min_depth, kmer_length, msa_kmers, cores):
	try:
		inputs = []
		tmp_fasta_output = output_prefix + '.fastq'
		for r in illumina_read_files:
			r = r.strip()
			if os.path.isfile(r):
				if is_gz_file(r):
					inputs.append(r)
					cat_cmd = "zcat %s > %s" % (r, tmp_fasta_output)
					subprocess.call(cat_cmd, shell=True)
				else:
					inputs.append(r)
					cat_cmd = "cat %s > %s" % (r, tmp_fasta_output)
					subprocess.call(cat_cmd, shell=True)
		if len(inputs) > 0:
			cmd1 = 'shopt -s extglob; jellyfish count -o %s -m %d -s 4000M -t %d -C %s' % (output_prefix + '.jf', kmer_length, cores, tmp_fasta_output)
			cmd2 = 'jellyfish dump -L %d -c -t %s -o %s' % (min_depth, output_prefix + '.jf', output_prefix + '.kmers.tab')

			subprocess.call(cmd1, shell=True)
			assert(os.path.isfile(output_prefix + '.jf'))
			subprocess.call('rm -f ' + tmp_fasta_output, shell=True)
			subprocess.call(cmd2, shell=True)
			assert(os.path.isfile(output_prefix + '.kmers.tab'))

			sample_kmers = set([])
			with open(output_prefix + '.kmers.tab') as oskf:
				for line in oskf:
					line = line.strip()
					ls = line.split('\t')
					if ls[0] in msa_kmers:
						sample_kmers.add(ls[0].upper())
						sample_kmers.add(revcomp(ls[0].upper()))
			return sample_kmers
	except:
		raise RuntimeError("Issues running JellyFish.")

def revcomp(seqx):
	seqxrc = Seq(seqx).reverse_complement()
	return (str(seqxrc))

def is_gz_file(filepath):
	with open(filepath, 'rb') as test_f:
		return binascii.hexlify(test_f.read(2)) == b'1f8b'

def is_fasta(fasta):
	try:
		with open(fasta) as of:
			SeqIO.parse(of, 'fasta')
		return True
	except:
		return False

def logParameters(parameter_names, parameter_values):
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		sys.stderr.write(pn + ': ' + str(pv) + '\n')

def logParametersToFile(parameter_file, parameter_names, parameter_values):
	parameter_handle = open(parameter_file, 'w')
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		parameter_handle.write(pn + ': ' + str(pv) + '\n')
	parameter_handle.close()

def logParametersToObject(logObject, parameter_names, parameter_values):
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		logObject.info(pn + ': ' + str(pv))