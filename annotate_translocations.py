# Rachel Kositsky
# 2019-04-16, 2019-04-22, 2019-04-24

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


def main(args):
	"""Produce annotated translocation table"""

	####################################################################
	### Read in translocation coordinates and initialize annotations ###
	####################################################################

	df = pd.read_csv(args.translocation_tsv, sep="\t", header=0)

	# Update these as you add more annotations
	# TODO: change this to be read from arguments instead of defined here
	annotation_types = ["Gene", "Promoter", "Panel_Category"]

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
		"-b", args.gene_bed, args.promoter_bed, args.panel_bed, "-names",
		"Gene", "Promoter", "Panel_Category"]

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
			a = BedtoolsAnnotation(line)
			# e.g. Gene_2
			col_name = a.annotation_category + "_" + a.location_order

			# Add on a 
			current_annotation = df.loc[a.df_row, col_name]
			if current_annotation:
				df.loc[a.df_row, col_name] = current_annotation + "," + a.annotation
			else:
				df.loc[a.df_row, col_name] = a.annotation

	# Delete temporary files here
	os.remove(annot_tmp_bed)
	os.remove(transloc_tmp_bed)

	######################################################
	### Rearrange annotations into tab-delimited table ###
	######################################################

	# Write out df as a tab-delimited csv
	df.to_csv(path_or_buf=args.out_tsv, sep="\t", index=False)


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
		help="BED file with category annotations (MYC, BCL2, BCL6,..)")

	# TODO: figure out some way to get a variable length list of pairs of names
	# and of BED file locations
	# e.g. -I (input) Gene genes.bed Promoter promoter.bed Enhancer enhancers.bed

	results = parser.parse_args(args)

	return results


if __name__ == '__main__':
    main(parse_args(sys.argv[1:]))