# TranslocationFinder

Scripts to find translocations using discordant reads from short-read sequencing aligned files (BAM files).


## Requirements

Python packages:
* pandas
* pybedtools
* pysam

Linux applications:
* libz-dev (for compiling pybedtools)
* bedtools
* samtools

## Running

1) Run select_discordant_reads.bash to get a BAM file and index with only discordant reads.

2) Run find_translocations.py to output a table with all translocations supported by the discordant reads provided.

3) Run annotate_translocations.py with annotation regions to annotate the translocations found.

Annotation regions are provided by BED files. Current annotations are gene-level, promoter-level, and category-level. 
Empty files can be provided if no annotation is desired for one of the categories.

## Authors

* [Rachel Kositsky](https://github.com/rkositsky)
