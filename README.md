# fSSU_pipeline

Run the pipeline on test data:

```
./fSSU_pipeline_v1.0.sh -t
```

Run the pipeline on test data with example files:

```
./fSSU_pipeline_v1.0.sh -1 testdata/mock_R1.fq -2 testdata/mock_R2.fq -b testdata/barcodes.txt -a CLC_assembly_cell -d testdata/mock16S.alignment -c 2
```

Run pipleline on real data:

```
./fSSU_pipeline_v1.0.sh -1 R1.fq -2 R2.fq -b barcodes.txt -a CLC_assembly_cell -d SILVA_128_94.alignment -c 2

where:
    -h  Show this help text
    -1  Read1 fastq file
    -2  Read2 fastq file
    -b  Sample barcodes
    -a  Path to CLC Assembly Cell installation
    -d  Path to aligned SSU reference database in fasta format with
        50000 alignment positions.
    -c  Number of CPUs to use
    -t  Flag to run pipeline on test data
```

Requirements - if your programs are named differently, then change the reference in the respective bash scripts:

```
mothur (version 1.37.6; https://www.mothur.org/wiki/Download_mothur)
parallel (version 20150422; https://www.gnu.org/software/parallel/)
CLC assembly cell (version 5.0.3; https://www.qiagenbioinformatics.com/products/clc-assembly-cell/)
perl (v5.18.2; https://www.perl.org/get.html)
```

# fSSU_pipeline_nanopore
Experimental pipeline for error-correction of Oxford Nanopore data using molecular tagging
