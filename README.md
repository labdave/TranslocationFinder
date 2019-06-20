# TranslocationFinder

Scripts to find translocations using discordant reads from short-read sequencing aligned files (BAM files).


## Requirements

Requires Python >3.5.

Python packages:
* pandas 0.24.2
* pybedtools 0.8.0
* pysam 0.15.2

Linux applications:
* libz-dev (for compiling pybedtools)
* bedtools
* samtools

## Running

Running Main.py runs the following scripts:

1) select_discordant_reads.bash: get a BAM file and index with only discordant reads
2) find_translocations.py: output a table with all translocations supported by the discordant reads provided
3) annotate_translocations.py: annotate translocations found with annotation regions

Annotation regions are provided by BED files. Current annotations are gene-level, promoter-level, category-level, and literature-derived.
Empty files can be provided if no annotation is desired for one of the categories.

## Authors

* [Rachel Kositsky](https://github.com/rkositsky)
