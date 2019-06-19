# Rachel Kositsky
# 2019-04-16, 2019-04-22, 2019-04-24, 2019-05-10, 2019-06-19

import argparse
import pandas as pd
import os
import subprocess
import sys


class BedtoolsAnnotation(object):
	def __init__(self, line):
		"""Parse a bedtools annotation line"""
		fields = line.strip("\n").split("\t")

		row, loc_order = fields[3].split("_")

		# The row number (0, 1, 2, ...)
		self.df_row = int(row)

		# The translocation number: "1" or "2"
		self.location_order = loc_order

		# e.g. Gene
		self.annotation_category = fields[4]

		# The value of the annotation, e.g. "MYC"
		self.annotation = fields[8]


def annotate_translocations(translocation_tsv, out_tsv, gene_bed, promoter_bed,
	panel_bed, literature_bed):
	"""Produce annotated translocation table

	Args:
		translocation_tsv: Input file with tab-separated translocation coordinates
		out_tsv: Output file path for annotated translocations
		gene_bed: BED file of gene body regions, annotated with gene names
		promoter_bed: BED file of promoter regions (2000 bp upstream of gene
			body), annotated with gene names
		panel_bed: BED file with category annotations (e.g. MYC-rearrangement
		literature_bed: BED file with literature-derived regions
	Returns: None
	Writes: table at out_tsv
	"""

	####################################################################
	### Read in translocation coordinates and initialize annotations ###
	####################################################################

	df = pd.read_csv(translocation_tsv, sep="\t", header=0)

	# Update these as you add more annotations
	# TODO: change this to be read from arguments instead of defined here
	annotation_types = ["Gene", "Promoter", "Panel_Category", "Literature"]

	# Add empty strings for all annotation types for both breakpoint ends
	# e.g. Gene_1, Gene_2, ...
	for a_t in annotation_types:
		for i in ["_1", "_2"]:
			col_name = a_t + i
			df[col_name] = ""

	##############################################
	### Run annotation with bedtools intersect ###
	##############################################

	# Make a temporary BED file for intersection with annotation files
	transloc_tmp_bed = "transloc.tmp.bed"
	with open(transloc_tmp_bed, "w") as out_f:
		# Separate out each breakpoint on its own line
		for i in range(df.shape[0]):
			# Write location 1 (0:3) and location 2 (3:6)
			# yes, it's mixing 0-based and 1-based counting; we need it for 
			# human readability later
			for j in range(1,3):
				coords = df.iloc[i, 3*(j-1):(3*(j-1)+3)].to_list()
				out_f.write("\t".join(map(str, coords)))
				out_f.write("\t{0}_{1}\n".format(i,j))


	# Intersect with all annotation BED files
	# TODO: update command to change with annotation arguments passed in
	cmd = ["bedtools", "intersect", "-wa", "-wb", "-a", transloc_tmp_bed, 
		"-b", gene_bed, promoter_bed, panel_bed,
		literature_bed, "-names"]
	# Tack on the annotation names to the command
	cmd = cmd + annotation_types

	annot_tmp_bed = "annot.tmp.bed"
	with open(annot_tmp_bed, "w") as out_f:
		annot_proc = subprocess.Popen(cmd, stdout=out_f)

		# Allow intersect_proc to receive a SIGPIPE if merge_proc exits
		#intersect_proc.stdout.close()
		# Run annot_proc
		annot_proc.communicate()

	# Read in annotations and add to dataframe
	# Note: 2+ annotation files required for the columns to work out
	with open(annot_tmp_bed, "r") as in_f:
		for line in in_f.readlines():
			ann = BedtoolsAnnotation(line)
			# e.g. Gene_2
			col_name = ann.annotation_category + "_" + ann.location_order

			# Add on annotation coulmn
			current_annotation = df.loc[ann.df_row, col_name]
			if current_annotation:
				df.loc[ann.df_row, col_name] = current_annotation + "," + ann.annotation
			else:
				df.loc[ann.df_row, col_name] = ann.annotation

	######################################################
	### Rearrange annotations into tab-delimited table ###
	######################################################

	# Move "Read" column last
	df = df.reindex(list([a for a in df.columns if a != 'Read_Pair_IDs']) + 
		['Read_Pair_IDs'], axis=1)

	# Sort dataframe to aid readability later
	df = df.sort_values(by=["N_Read_Pairs", "Chrom_1", "Start_1", "Chrom_2", 
		"Start_2"], ascending=[False, True, True, True, True])

	# Write out df as a tab-delimited csv
	df.to_csv(path_or_buf=out_tsv, sep="\t", index=False)

	# Delete temporary files
	os.remove(annot_tmp_bed)
	os.remove(transloc_tmp_bed)


def parse_args(args=None):
	"""Parse command line arguments"""
	parser = argparse.ArgumentParser(
		description="Annotate translocations from find_translocations.py")

	parser.add_argument("translocation_tsv",
		help="Input file with tab-separated translocation coordinates")

	parser.add_argument("out_tsv",
		help="Output file path for annotated translocations")

	parser.add_argument("gene_bed", 
		help="BED file of gene body regions, annotated with gene names")

	parser.add_argument("promoter_bed", 
		help="BED file of promoter regions (2000 bp upstream of gene body), annotated with gene names")

	parser.add_argument("panel_bed", 
		help="BED file with category annotations (e.g. MYC-rearrangement)")

	parser.add_argument("literature_bed",
		help="BED file with literature-derived regions")

	# TODO: figure out some way to get a variable length list of pairs of names
	# and of BED file locations
	# e.g. -I (input) Gene genes.bed Promoter promoter.bed Enhancer enhancers.bed

	results = parser.parse_args(args)

	return results


if __name__ == '__main__':
	args = parse_args(sys.argv[1:])
	annotate_translocations(args.translocation_tsv, args.out_tsv, args.gene_bed,
		args.promoter_bed, args.panel_bed, args.literature_bed)
