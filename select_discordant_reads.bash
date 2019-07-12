#!/bin/bash

# Rachel Kositsky
# Created: 2019-04-11
# Updated: 2019-06-20

# Goal: given a BAM file, produce another BAM file that selects discordant reads

# Usage if don't have inputs
if [[ $# -ne 4 ]] ; then
    echo "Usage: select_discordant_reads.bash [in.bam] [out.bam] [nr_cpus] [min_map_quality]"
    echo "Given a BAM file, produce another BAM file that selects discordant reads."
    echo "  in.bam, out.bam: input and output file paths"
    echo "  nr_cpus: number of CPUs over which to parallelize"
    echo "  min_map_quality: minimum mapping quality for selected reads"
    exit 0
fi

input=$1
output=$2
nr_cpus=$3
min_mq=$4

tmp_output="tmp.bam"

# 1) Extract discordant reads
samtools view -b -@ ${nr_cpus} -f 1 -F 1038 -q ${min_mq} -o ${tmp_output} ${input}
samtools index ${tmp_output}

# 2) Remove duplicates. -S: consider paired-end reads as single-end.
samtools rmdup -S ${tmp_output} ${output}
samtools index ${output}
