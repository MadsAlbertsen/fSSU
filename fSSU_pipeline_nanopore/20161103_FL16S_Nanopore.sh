#############################################################################
# 									                               #
# Shell script for generating error corrected FL16S and getting error rates #
#									                               #
# Use at your own RISK!							                     #
#############################################################################

####################
#     Variables    #
####################
ID_adapt=0.1;
ID_cluster=0.8;
LINKtoCANU=/space/users/rkirke08/Desktop/canu/canu-1.3/Linux-amd64/bin;
# Update path to include poretools installation
export PATH=$PATH:/space/users/rkirke08/.local/bin

####################
# End of variables #
####################

###############################################
# Depends on the following files and software #
###############################################
# folder with fast5 files "data/pass/"
# file with reference 16S sequences "mockrRNAall.fasta"
# perl script "F16S.cluster.split.pl"
# poretools
# cutadapt
# usearch8.1
# CANU
###############################################
# End of dependencies			      #
###############################################

# Extract fastq files
poretools fastq --type 2D data/pass/ > data/2D.fq

# Rename headers (Some tools do not accept the long poretools headers)
awk '{print (NR%4 == 1) ? "@" ++i : $0}' data/2D.fq | sed -n '1~4s/^@/>/p;2~4p' > 2Dr.fa

# Find adapters
cutadapt -g AAAGATGAAGAT -e $ID_adapt -O 12 -m 1300 --untrimmed-output un1.fa -o a1.fa 2Dr.fa
cutadapt -a ATGGATGAGTCT -e $ID_adapt -O 12 -m 1300 --discard-untrimmed -o a1_a2.fa a1.fa

usearch8.1 -fastx_revcomp un1.fa -label_suffix _RC -fastaout un1_rc.fa
cutadapt -g AAAGATGAAGAT -e $ID_adapt -O 12 -m 1300 --discard-untrimmed -o ua1.fa un1_rc.fa
cutadapt -a ATGGATGAGTCT -e $ID_adapt -O 12 -m 1300 --discard-untrimmed -o ua1_a2.fa ua1.fa

cat a1_a2.fa ua1_a2.fa > c.fa

# Extract barcodes
cut -c1-12 c.fa > i1.fa
rev c.fa | cut -c1-12 | rev > i2.fa

paste i1.fa i2.fa -d "" | cut -f1-2 -d ">" > i1i2.fa

# Cluster barcodes
usearch8.1 -cluster_fast i1i2.fa -id $ID_cluster -centroids nr.fa -uc res.uc -sizeout

# Extract raw sequences
perl F16S.cluster.split.pl -c res.uc -i c.fa -m 3 -f 50 -r 40
# Count number of files in directory
find clusters -type f | wc -l

FILES=clusters/*.fa
for OTU in $FILES

do
	wc -l $OTU >> lines.txt	
done

FILES=clusters/*.fa
for OTU in $FILES

do
	OTUNO=$(echo $OTU | cut -f2 -d\/);
	# Rename header
	sed "s/>/>$OTUNO/" clusters/$OTUNO > clusters/newHeaders_$OTUNO

	# Correct reads using CANU
	$LINKtoCANU/canu -correct -p OTU_$OTUNO -d clusters/OTU_$OTUNO genomeSize=1.5k -nanopore-raw  $OTU

	# Unsip corrected reads
	gunzip clusters/OTU_$OTUNO/OTU_$OTUNO.correctedReads.fasta.gz

	sed -i "s/>/>$OTUNO/" clusters/OTU_$OTUNO/OTU_$OTUNO.correctedReads.fasta

	# Call consensus using Usearch
	usearch8.1 -cluster_fast clusters/OTU_$OTUNO/OTU_$OTUNO.correctedReads.fasta -id 0.9 -centroids clusters/OTU_$OTUNO/nr_cor_$OTUNO.fa -uc clusters/OTU_$OTUNO/res_$OTUNO.uc -sizeout -consout clusters/OTU_$OTUNO/Ucons_CANUcor_$OTUNO.fa

	sed -i "s/>/>$OTUNO/" clusters/OTU_$OTUNO/Ucons_CANUcor_$OTUNO.fa

	# Map reads to references to estimate error rate
	# Raw reads
	# Map FL16S back to references
	usearch8.1 -usearch_global clusters/newHeaders_$OTUNO -db mockrRNAall.fasta -strand both -id 0.60 -top_hit_only -maxaccepts 10 -query_cov 0.5 -userout clusters/map_raw_$OTUNO.txt -userfields query+target+id+ql+tl+alnlen

	# Usearch consensus corrected sequence
	# Map FL16S back to references
	usearch8.1 -usearch_global clusters/OTU_$OTUNO/Ucons_CANUcor_$OTUNO.fa -db mockrRNAall.fasta -strand both -id 0.60 -top_hit_only -maxaccepts 10 -query_cov 0.5 -userout clusters/map_cor_Ucons_$OTUNO.txt -userfields query+target+id+ql+tl+alnlen

	cat clusters/map_raw_$OTUNO.txt >> myfile.txt
	cat clusters/map_cor_Ucons_$OTUNO.txt >> myfile.txt

	# Collect corrected sequences
	cat clusters/OTU_$OTUNO/Ucons_CANUcor_$OTUNO.fa >> final_corrected.fa
done
