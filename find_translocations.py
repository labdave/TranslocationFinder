# Identify translocations given discordant reads, optionally filtering by a BED
# with target regions
#
# Rachel Kositsky
# Created: 2019-04-14

import argparse
import os
import pybedtools
import pysam
import subprocess
import sys
import time


class Translocation(object):
	"""Class for translocations"""
	def __init__(self, read_1, read_2):
		"""Initializes Translocation object
		Args:
			read_1: pysam AlignmentSegment with the earlier chromosome; not
					necessarily read1
			read_2: pysam AlignmentSegment with the later chromosome; not
					necessarily read2
		"""

		# Location 1 is set from initial read that lends its support;
		# expand when another supporting read is nearby
		self.chrom_1 = read_1.reference_name
		self.start_1 = read_1.reference_start
		self.end_1 = read_1.reference_start + read_1.reference_length

		# Location 2 is not set yet; initialize with the first paired read,
		# but be ready to expand it if another supporting read nearby is found
		self.chrom_2 = read_2.reference_name
		self.start_2 = read_2.reference_start
		self.end_2 = read_2.reference_start + read_2.reference_length

		# Just keep track of read IDs, not the reads themselves
		self.read_ids = [read_1.qname]
		self.n_read_pairs = 1

	def add_read_pair(self, read_1, read_2):
		"""Adds a supporting read pair to the translocation
		Args:
			read_1: pysam AlignmentSegment with the earlier chromosome; not
					necessarily read1
			read_2: pysam AlignmentSegment with the later chromosome; not
					necessarily read2
		"""

		# Chromosome should be the same
		assert(read_1.reference_name == self.chrom_1)
		assert(read_2.reference_name == self.chrom_2)

		# Update start and ending positions
		self.start_1 = min(self.start_1, read_1.reference_start)
		self.end_1 = max(self.end_1, read_1.reference_start + read_1.reference_length)

		self.start_2 = min(self.start_2, read_2.reference_start)
		self.end_2 = max(self.end_2, read_2.reference_start + read_2.reference_length)

		# Update read IDs and number of supporting read pairs
		self.read_ids.append(read_1.qname)
		self.n_read_pairs += 1
	
	def near(self, read, max_distance):
		"""Returns True/False based on whether read is within max_distance
		of the endpoint of the translocation"""
		
		if (read.reference_name != self.chrom_2):
			return False

		start = read.pos
		end = read.pos + read.rlen

		result = ((abs(self.start_2 - start) <= max_distance) and 
			(abs(self.end_2 - end) <= max_distance))

		return result

	def get_output_line(self):
		"""Output line for translocation file."""

		out_list = [self.chrom_1, self.start_1, self.end_1, self.chrom_2, 
			self.start_2, self.end_2, self.n_read_pairs, 
			",".join(self.read_ids)]
		out_list = map(str, out_list)

		return "\t".join(out_list) + "\n"


class MergedLocation(object):
	"""Class for intermediate result of merged lines"""

	def __init__(self, line):
		"""Input: Output from bedtools merge call"""
		fields = line.strip("\n").split("\t")
		self.chrom = fields[0]
		self.start = int(fields[1])
		self.end = int(fields[2])
		self.n_reads = int(fields[3])
		self.read_ids = fields[4]

	def __str__(self):
		return "{0}:{1}-{2}".format(self.chrom, self.start, self.end)


def setup(out_dir, discordant_reads_bam):
	"""Set up for find_translocations

	Args:
		out_dir: output directory
		discordant_reads_bam: input BAM file as a pysam object
	Returns:
		start_time: the starting time of the program
		ref_order: a list of chromosomes in reference order
		ref_key: function to get the order of a given chromosome
	"""
	start_time = time.time()

	# Create output directory if it does not exist
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)

	# Parse in reference contigs (chromosome) ordering from input BAM file
	header_dict = discordant_reads_bam.header.as_dict()
	ref_order = []

	for sequence in header_dict["SQ"]:
		ref_order.append(sequence["SN"])

	ref_key = lambda ref: ref_order.index(ref)

	return (start_time, ref_order, ref_key)


def find_translocations(discordant_reads_bam, out_dir,
	BED_filter=None, merge_distance=1000, min_read_pairs=2,
	min_mapping_quality=30):
	"""Identify translocations given discordant reads,
	optionally filtering by a BED with target region of one end.

	How: use bedtools merge to get pile ups of reads in one location,
	then filter through the pile ups to see if you get any with pile ups in 
	the second location.

	Args:
	  discordant_reads_bam: BAM file with discordant reads
	  out_dir: Output directory
	  BED_filter: BED file on which to filter
	  merge_distance: Maximum distance between reads for which they're still
	    considered part of the same translocation
	  min_read_pairs: Minimum supporting reads for a translocation to be called
	  min_mapping_quality: Minimum mapping quality for any supporting reads
	Returns: None
	Writes in output directory:
	  merged_reads.bed: intermedite file with merged discordant reads
	  translocations.tsv: output file with translocations
	"""

	###########
	## Setup ##
	###########

	bamfile = pysam.AlignmentFile(discordant_reads_bam, "rb")
	(start_time, ref_order, ref_key) = setup(out_dir, bamfile)

	###############################################
	## Pile up reads to get first breakpoint end ##
	###############################################

	# Use bedtools merge
	# Issue: can't merge the mates along with the original guys.
	# Workaround: just use this for merging R1s, then filter later on R2s
	merged_file = os.path.join(out_dir, "merged_reads.bed")

	# Convert BAM to BED
	bamtobed_proc = subprocess.Popen(["bedtools", "bamtobed", 
		"-i", discordant_reads_bam], stdout=subprocess.PIPE)

	# If a BED target file was specified, filter the reads before merging
	if BED_filter:
		# Filter using intersect
		# run with -nonamecheck to avoid a warning from GL000216.2 chromosome
		intersect_proc = subprocess.Popen(["bedtools", "intersect", "-a", "-", 
				"-b", BED_filter, "-nonamecheck"],
				stdin=bamtobed_proc.stdout, stdout=subprocess.PIPE)

		# Got error from bedtools merge about unsorted input, so sort here
		sort_proc = subprocess.Popen(["bedtools", "sort", "-i", "-"],
			stdin=intersect_proc.stdout, stdout=subprocess.PIPE)

		# Merge the reads within leeway distance
		with open(merged_file, "w") as out_f:
			merge_proc = subprocess.Popen(["bedtools", "merge", "-i", "-", 
				"-d", str(merge_distance), "-c", "4", "-o", "count,distinct"],
				stdin=sort_proc.stdout, stdout=out_f)

			# Allow bamtobed_proc to receive a SIGPIPE if intersect_proc exits
			bamtobed_proc.stdout.close()
			# Allow intersect_proc to receive a SIGPIPE if sort_proc exits
			intersect_proc.stdout.close()
			# Allow sort_proc to receive a SIGPIPE if merge_proc exits
			sort_proc.stdout.close()
			# Run merge_proc
			merge_proc.communicate()

		# Set up pybedtools for later
		bed_filter_file = pybedtools.BedTool(BED_filter)

	# Otherwise just merge reads without the filtering step beforehand
	else:
		# Got error from bedtools merge about unsorted input, so sort here
		sort_proc = subprocess.Popen(["bedtools", "sort", "-i", "-"],
			stdin=bamtobed_proc.stdout, stdout=subprocess.PIPE)

		# Merge the reads within leeway distance
		merged_file = os.path.join(out_dir, "merged_reads.bed")
		with open(merged_file, "w") as out_f:
			merge_proc = subprocess.Popen(["bedtools", "merge", "-i", "-", 
				"-d", str(merge_distance), "-c", "4",
				"-o", "count,distinct"], 
				stdin=sort_proc.stdout, stdout=out_f)

			# Allow bamtobed_proc to receive a SIGPIPE if sort_proc exits
			bamtobed_proc.stdout.close()
			# Allow sort_proc to receive a SIGPIPE if merge_proc exits
			sort_proc.stdout.close()
			# Run merge_proc
			merge_proc.communicate()

	#############################################
	## Go through merged lines and compare R2s ##
	#############################################

	bamfile = pysam.AlignmentFile(discordant_reads_bam, "rb")
	tranlocation_file = os.path.join(out_dir, "translocations.tsv")
	n_merged = 0
	n_written = 0

	with open(merged_file, "r") as in_f, open(tranlocation_file, "w") as out_f:
		# Write header for translocation output
		out_f.write("\t".join(["Chrom_1", "Start_1", "End_1",
			"Chrom_2", "Start_2", "End_2", "N_Read_Pairs", "Read_Pair_IDs"]))
		out_f.write("\n")

		# Loop through merged lines that are possible translocation locations
		for line in in_f.readlines():
			n_merged += 1

			# The first location of a possible translocation
			loc1 = MergedLocation(line)

			# Skip over if not enough read support in first location
			if loc1.n_reads < min_read_pairs:
				continue

			# List of all translocations for this potential translocation start
			translocation_list = []

			## Get locations of mate pairs ##
			# TODO: include chimeric reads/supplementary alignments
			n_mates = 0
			for read in bamfile.fetch(loc1.chrom, loc1.start, loc1.end):
				# Skipping over check here for read ID being part of region

				# Keep track of number reads in that region lining up
				n_mates += 1

				# Get paired read of the one in the original pileup
				try:
					mate = bamfile.mate(read)
				except ValueError:
					continue

				# Filter on mapping quality
				if mate.mapq < min_mapping_quality: continue

				mate_ref = mate.reference_name
				read_ref = read.reference_name

				# Only consider translocations with known chromosomes
				if mate_ref not in ref_order: continue

				# If BED file not provided, only consider translocations for 
				# chromosomes ahead of you
				if BED_filter == None:
					if ref_key(mate_ref) <= ref_key(read_ref): continue
				# If BED file provided, do inter-chromosomal for everybody, but
				# for BED file mates, only consider translocations for
				# chromosomes ahead of you so that you don't double count
				# translocations
				else:
					# We only capture inter-chromosomal breakpoints here; 
					# ignore any read pairs on the same chromosome
					if mate_ref == read_ref:
						continue
					else:
						# If it's in the BED target file, go back to WGS rules;
						# only keep if it's in a larger chromosome to not 
						# double-print
						mate_interval = pybedtools.Interval(mate_ref,
							mate.reference_start, 
							mate.reference_start + mate.reference_length)
						# .any_hits is 1 when there's an overlap
						if bed_filter_file.any_hits(mate_interval):
							if ref_key(mate_ref) <= ref_key(read_ref): continue

				## Add a translocation ##
				# If we got this far, we want to use this read pair
				# Try to add the mate to an existing translocation
				for t in translocation_list:
					if t.near(mate, merge_distance):
						t.add_read_pair(read, mate)
						break
				# If no translocation fits, then make a new one
				else:
					t = Translocation(read, mate)
					translocation_list.append(t)

			# Check that number of reads you trawled for mates matches bedtools
			# output
			if n_mates != loc1.n_reads:
				raise Exception("Reads fetched in region different from merged "
					"number of reads: merged {0}, but fetched {1} in region "
					"{2}".format(loc1.n_reads, n_mates, loc1))
				
			# Sort translocation list by chromosome, then by starting position
			translocation_list = sorted(translocation_list, 
				key=lambda t: t.start_2)
			translocation_list = sorted(translocation_list, 
				key=lambda t: ref_key(t.chrom_2))

			# Write any translocations which have enough reads
			for t in translocation_list:
				if t.n_read_pairs >= min_read_pairs:
					out_f.write(t.get_output_line())
					n_written += 1

	bamfile.close()
	print("----------------------")
	print("Processed {0} merged lines".format(n_merged))
	print("Found {0} translocations, and wrote them at {1}".format(n_written,
		tranlocation_file))

	finish_time = time.time() - start_time
	print("Time: {0:.0f} sec / {1:.0f} min / {2:.0f} hours".format(
		finish_time, finish_time/60., finish_time/3600.))


def parse_args(args=None):
	"""Parse command line arguments"""
	parser = argparse.ArgumentParser(
		description="Identify translocations given discordant reads")

	parser.add_argument("discordant_reads_bam",
		help="BAM file with discordant reads")

	parser.add_argument("out_dir",
		help="Output directory")

	parser.add_argument("-B", "--BED_filter", 
		help="BED file on which to filter", default=None)

	parser.add_argument("-D", "--merge_distance", default=1000, type=int,
		help="Maximum distance between reads for which they're still considered"
		" part of the same translocation (default 1000 bp)")

	parser.add_argument("-R", "--min_read_pairs", default=2, type=int,
		help="Minimum supporting reads for a translocation to be called"
		" (default 2 read pairs)")

	parser.add_argument("-Q", "--min_mapping_quality", default=30, type=int,
		help="Minimum mapping quality for any supporting reads (default 30)")

	results = parser.parse_args(args)

	return results


if __name__ == '__main__':
	a = parse_args(sys.argv[1:])
	find_translocations(a.discordant_reads_bam, a.out_dir, a.chromosome_list,
		a.BED_filter, a.merge_distance, a.min_read_pairs, a.min_mapping_quality)
