#!/usr/local/bin/python3

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
		gene_bed: BED file of gene body regions, annotated with gene names
		promoter_bed: BED file of promoter regions (2000 bp upstream of gene
			body), annotated with gene names
		target_bed: BED file with annotations from target panel (e.g. MYC-exon1)
		literature_bed: BED file with literature-derived regions
	"""

	# If BED file given, use that, otherwise use literature BED file as default.
	# -1 means no filtering.
	r_dir = os.path.join(dir_path, "resources")

	if args.BED_filter:
		if args.BED_filter == "-1":
			BED_filter = None
		else:
			BED_filter = args.BED_filter
	else:
		BED_filter = os.path.join(r_dir, "both_panels_MYC_BCL2_BCL6.bed")

	gene_bed = args.gene_bed if args.gene_bed else os.path.join(r_dir, "genes.bed")
	promoter_bed = args.promoter_bed if args.promoter_bed else os.path.join(r_dir, "promoters.bed")
	target_bed = args.target_bed if args.target_bed else os.path.join(r_dir, "panel.gencode.bed")
	literature_bed = args.literature_bed if args.literature_bed else os.path.join(r_dir, "literature.bed")

	return (BED_filter, gene_bed, promoter_bed, target_bed, literature_bed)


def main(args):
	# Get folder of current executed file to get the repository path
	dir_path = os.path.dirname(os.path.realpath(__file__))

	# Create temporary output directory
	out_dir = "translocation_work_dir"
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)

	# Parse arguments and fill in defaults if needed
	(BED_filter, gene_bed, promoter_bed,
		target_bed, literature_bed) = define_defaults(dir_path, args)

	# Call select_discordant_reads. Its output BAM file will be written in the
	# output directory.
	print("Extracting discordant reads...")
	extraction_script = os.path.join(dir_path, "select_discordant_reads.bash")
	out_bam = os.path.join(out_dir, "discordant_reads.bam")
	subprocess.call([extraction_script, args.in_bam, out_bam])

	# Call find_translocations.
	print("Finding translocations...")
	find_translocations.find_translocations(out_bam, out_dir, 
		BED_filter, args.merge_distance, args.min_read_pairs, 
		args.min_mapping_quality)

	# Call annotation file. Use files from resources/ defined earlier.
	print("Annotating translocations...")
	translocation_tsv = os.path.join(out_dir, "translocations.tsv")
	annotate_translocations.annotate_translocations(translocation_tsv,
		args.out_path, gene_bed, promoter_bed, target_bed, literature_bed)

	print("Done!")


def parse_args(args=None):
	"""Parse command line arguments"""
	parser = argparse.ArgumentParser(
		description="Find and annotate translocations using discordant reads "
		"from short-read sequencing aligned files")

	parser.add_argument("in_bam", help="Input BAM file")

	parser.add_argument("out_path",
		help="Path to output table with annotated translocations")

	parser.add_argument("--BED_filter", type=str,
		help="BED file where translocations can be called, or -1 for no "
		"filtering. Default: MYC, BCL2, BCL6, and other genes.",
		default=None)

	parser.add_argument("--merge_distance", default=1000, type=int,
		help="Maximum distance between reads where they're still considered"
		" part of the same translocation. Default: 1000 bp")

	parser.add_argument("--min_read_pairs", default=2, type=int,
		help="Minimum supporting reads for a translocation to be called"
		" . Default: 2 read pairs")

	parser.add_argument("--min_mapping_quality", default=30, type=int,
		help="Minimum mapping quality for supporting reads. Default: MQ >= 30")

	parser.add_argument("--gene_bed",
		help="BED file of gene body regions, annotated with gene names")

	parser.add_argument("--promoter_bed",
		help="BED file of promoter regions (2000 bp upstream of gene body), "
		"annotated with gene names")

	parser.add_argument("--target_bed",
		help="BED file with annotations from target panel (e.g. MYC-exon1)")

	parser.add_argument("--literature_bed",
		help="BED file with literature-derived regions")

	results = parser.parse_args(args)

	return results

if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))
