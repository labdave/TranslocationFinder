#!/bin/bash

# Rachel Kositsky
# Created: 2019-04-11
# Updated: 2019-04-11

# Goal: given a BAM file, produce another BAM file that selects discordant reads

# Usage if don't have inputs
if [[ $# -ne 2 ]] ; then
    echo "Usage: select_discordant_reads.bash [in.bam] [out.bam]"
    echo "Given a BAM file, produce another BAM file that selects discordant reads"
    exit 0
fi

input=$1
output=$2

 samtools view -b -@ 4 -f 1 -F 1038 -q 30 -o ${output} ${input}

 samtools index ${output}
