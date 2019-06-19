# Get discordant reads, identify translocations, and annotate.
#
# Rachel Kositsky
# Created: 2019-06-17
# Updated: 2019-06-19

import annotate_translocations
import argparse
import find_translocations
import os
import subprocess
import sys

def define_defaults(dir_path, args):
	"""Defines default paths to resources when running main()
	
	Args:
		dir_path: directory of currently executed file/repository path
		args: command line parsed arguments
	Returns:
		chromosome_list: sorted chromosome names (GENCODE default)
		BED_filter: BED genomic region file for where to call translocations
	"""
	
	# # Use GENCODE chromosome names by default.
	# if args.chromosome_list:
	# 	chromosome_list = args.chromosome_list
	# else:
	# 	chromosome_list = os.path.join(dir_path, "resources", 
	# 		"GRCh38_gencode.txt")

	# If BED file given, use that, otherwise use literature BED file as default. 
	if args.BED_filter:
		BED_filter = args.BED_filter
	else:
		BED_filter = os.path.join(dir_path, "resources", 
			"both_panels_MYC_BCL2_BCL6.bed")

	return BED_filter


def main(args):
	# Get folder of current executed file to get the repository path
	dir_path = os.path.dirname(os.path.realpath(__file__))

	# Create temporary output directory
	out_dir = "translocation_work_dir"
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)

	# Call select_discordant_reads. Its output BAM file will be written in the
	# output directory.
	extraction_script = os.path.join(dir_path, "select_discordant_reads.bash")
	out_bam = os.path.join(out_dir, "discordant_reads.bam")
	subprocess.call([extraction_script, args.in_bam, out_bam])

	# Call find_translocations. 
	BED_filter = define_defaults(dir_path, args)
	find_translocations.find_translocations(out_bam, out_dir, 
		BED_filter, args.merge_distance, args.min_read_pairs, 
		args.min_mapping_quality)

	# Call annotation file. Use files from resources/.
	translocation_tsv = os.path.join(out_dir, "translocations.tsv")
	gene_bed = os.path.join(dir_path, "resources", "genes.bed")
	promoter_bed = os.path.join(dir_path, "resources", "promoters.bed")
	panel_bed = os.path.join(dir_path, "resources", "panel.gencode.bed")
	literature_bed = os.path.join(dir_path, "resources", "literature.bed")

	annotate_translocations.annotate_translocations(translocation_tsv, 
		args.out_path, gene_bed, promoter_bed, panel_bed, literature_bed)


def parse_args(args=None):
	"""Parse command line arguments"""
	parser = argparse.ArgumentParser(
		description="Find and annotate translocations using discordant reads "
		"from short-read sequencing aligned files")

	parser.add_argument("in_bam", help="Input BAM file")

	parser.add_argument("out_path",
		help="Path to output table with annotated translocations")

	parser.add_argument("-B", "--BED_filter",
		help="BED file where translocations can be called. "
		"Default: MYC, BCL2, BCL6, and other genes.",
		default=None)

	parser.add_argument("-D", "--merge_distance", default=1000, type=int,
		help="Maximum distance between reads where they're still considered"
		" part of the same translocation. Default: 1000 bp")

	parser.add_argument("-R", "--min_read_pairs", default=2, type=int,
		help="Minimum supporting reads for a translocation to be called"
		" . Default: 2 read pairs")

	parser.add_argument("-Q", "--min_mapping_quality", default=30, type=int,
		help="Minimum mapping quality for supporting reads. Default: MQ >= 30")

	results = parser.parse_args(args)

	return results

if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
