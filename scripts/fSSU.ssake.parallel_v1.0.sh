#!/bin/bash

# Parallel assembly of read bins

# Arguments
datafolder=$1
n_reads=${2:-1000}
cov_min=${3:-10}
overlap=${4:-16}
minbase=${5:-2}
baseratio=${6:-0.51}
outfolder=${7:-"assembly"}
contig_min=${8:-700}
error=${9:-20}
read_min=${10:-50}
a5=${11:-"AAAAAAAAAAAAA"}
a3=${12:-"GGGCAATATCAGCAC"}

# Preparation
t1=$(date +%s.%N)
mkdir $outfolder
args="$n_reads $cov_min $overlap $minbase $baseratio $outfolder $contig_min $error $read_min $a5 $a3"

### Assembly ###
find  $datafolder/ -name '*.fastq' |\
sed -e "s/$/ $args/" |\
parallel --progress -j90% 'sh scripts/fSSU.ssake_v1.0.sh'

### Stats ###
t2=$(date +%s.%N)
td=$(echo "$t2 - $t1" | bc)

echo "n_reads cov_min overlap minbase baseratio outfolder contig_min error read_min adp5 adp3 Process_Time" >> $outfolder/${outfolder}_stats.txt
echo "$args $td" >> $outfolder/${outfolder}_stats.txt
echo "## Assembly Stats ##" >> $outfolder/${outfolder}_stats.txt
for sf in $outfolder/*.stats
do
	cat $sf >> $outfolder/${outfolder}_stats.txt
done

### Extract longest contig from each assembly ###

for file in $outfolder/*.fa
do
  tag=$(basename $file | sed 's/\.fa//')
  head -n2 $file | sed "1s/^.contig./\>$tag/" >> $outfolder/all_sequences.fasta
done


