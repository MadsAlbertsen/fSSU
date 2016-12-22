#!/bin/bash
# Requirements - if your programs are named differently, then change the reference in the respective bash scripts
# cutadapt (version 1.10; http://cutadapt.readthedocs.io/en/stable/)
# barrnap (version 0.7; https://github.com/tseemann/barrnap)
# Usearch8.1 (version v8.1.1861; http://www.drive5.com/usearch/)
# SSAKE (version 3.8.4; https://github.com/warrenlr/SSAKE)

# Define input files and arguments
raw_r1=${1:-R1.fq}
raw_r2=${2:-R2.fq}

# Extract read-bins based on tags
perl scripts/fSSU.extract_v4.2.pl -f $raw_r1 -r $raw_r2 -q
mkdir readbins
mv -t readbins/ AnchorReads A1_Read_Stats.txt A2_Read_Stats.txt Anchor_Stats.txt stats.txt

# Assemble full-length sequences
bash scripts/fSSU.ssake.parallel_v1.0.sh readbins/AnchorReads/CTAGTACGATAGAGAG
mkdir results
cp assembly/all_sequences.fasta results/

# Filter, split (16S, 18S) and cluster (uses results/all_sequences.fasta)
bash scripts/fSSU.post.assembly_v1.0.sh
