#!/bin/bash

# Preparation
data=$(echo $1| cut -f1 -d" ")
n_reads=$(echo $1 |cut -f2 -d" ")
cov_min=$(echo $1 |cut -f3 -d" ")
overlap=$(echo $1 |cut -f4 -d" ")
minbase=$(echo $1 |cut -f5 -d" ")
baseratio=$(echo $1 |cut -f6 -d" ")
outfolder=$(echo $1 |cut -f7 -d" ")
contig_min=$(echo $1 |cut -f8 -d" ")
error=$(echo $1 |cut -f9 -d" ")
read_min=$(echo $1 |cut -f10 -d" ")
a5=$(echo $1 |cut -f11 -d" ")
a3=$(echo $1 |cut -f12 -d" ")

a5rc=$(echo $a5 | rev | tr ATGC TACG)
a3rc=$(echo $a3 | rev | tr ATGC TACG)
filename=$(basename $data | sed 's/\.fast.//')
mkdir $outfolder/$filename
workdir=${outfolder}/${filename}

### Assembly ###
# Subsample, trim and convert to fasta
head -n $(($n_reads * 4)) $data |\
cutadapt -a $a5 -a $a5rc -a $a3 -a $a3rc -q $error -m $read_min - --quiet |\
sed -n '1~4s/^@/>/p;2~4p' > $workdir/reads.fa

# SSAKE
SSAKE -f $workdir/reads.fa -w $cov_min -m $overlap -o $minbase -r $baseratio -p 0 -b $workdir/$filename -h -v default >/dev/null

# Stats
bp_max=$(cat $workdir/$filename.contigs | sed -n 'p;n' | cut -f2 -d"|" | tr -d size | sort -nr | head -n1)
contig_bp=$(cat $workdir/$filename.contigs | sed -n 'p;n' | cut -f2 -d"|" | tr -d size)
contig_n=0
for int in $contig_bp
do
  [ $int -gt 500 ]  && contig_n=$((contig_n+1))
done
echo "$filename $bp_max $contig_n" >> $outfolder/$filename.stats

### Cleanup ###
mv $workdir/$filename.contigs $outfolder/"$filename.fa"
#mv $workdir/cov.txt $outfolder/"$filename.cov"
rm -rf $workdir