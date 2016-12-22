Assembly=results/all_sequences.fasta
OutF=results

echo "Post assembly: Trimming"
cutadapt -g TGGTGCTGATATTGC -a GCAATATCAGCACCA -o $OutF/trimmed.fa -m 1200 $Assembly --quiet

echo "Post assembly: Classifying (Kingdom + rRNA type)"
barrnap $OutF/trimmed.fa --reject 0.3 --threads 40 --kingdom bac --quiet > $OutF/Bac.txt
barrnap $OutF/trimmed.fa --reject 0.3 --threads 40 --kingdom arc --quiet > $OutF/Arc.txt
barrnap $OutF/trimmed.fa --reject 0.3 --threads 40 --kingdom euk --quiet > $OutF/Euk.txt
sed -i 's/Name=/Name=B/g' $OutF/Bac.txt
sed -i 's/Name=/Name=A/g' $OutF/Arc.txt
cat $OutF/Bac.txt $OutF/Arc.txt $OutF/Euk.txt > $OutF/classified.txt
perl scripts/F16S.rRNA.split.pl -i $OutF/trimmed.fa -o $OutF -c $OutF/classified.txt

echo "Post assembly: Clustering and chimeral removal (Bacteria 16S rRNA)"
usearch8.1 -cluster_fast $OutF/B16S_rRNA.fa -id 0.97 -centroids $OutF/B16S_rRNA_cluster97.fa -sizeout -sort size -quiet
usearch8.1 -uchime_denovo $OutF/B16S_rRNA_cluster97.fa -nonchimeras $OutF/B16S_rRNA_nonchimeras_denovo97.fa -quiet

echo "Post assembly: Clustering and chimeral removal (Archaea 16S rRNA)"
usearch8.1 -cluster_fast $OutF/A16S_rRNA.fa -id 0.97 -centroids $OutF/A16S_rRNA_cluster97.fa -sizeout -sort size -quiet
usearch8.1 -uchime_denovo $OutF/A16S_rRNA_cluster97.fa -nonchimeras $OutF/A16S_rRNA_nonchimeras_denovo97.fa -quiet

echo "Post assembly: Clustering and chimeral removal (Eukaryot 18S rRNA)"
usearch8.1 -cluster_fast $OutF/18S_rRNA.fa -id 0.97 -centroids $OutF/18S_rRNA_cluster97.fa -sizeout -sort size -quiet
usearch8.1 -uchime_denovo $OutF/18S_rRNA_cluster97.fa -nonchimeras $OutF/18S_rRNA_nonchimeras_denovo97.fa -quiet

echo "Post assembly: Removing temp files"
rm $OutF/classified.txt $OutF/Bac.txt $OutF/Euk.txt $OutF/Arc.txt

echo "Done. Enjoy."