# TaG-EM Code and Files
This repository contains code and reference files related to the following publication:
Mendana, Jorge Blanco, Donovan, Margaret, O'Brien, Lindsey Gengelbach, Auch, Benjamin, Garbe, John, Gohl, Daryl M. (2025) Deterministic Genetic Barcoding for Multiplexed Behavioral and Single-Cell Transcriptomic Studies. eLife2025;12:RP88334 DOI: https://doi.org/10.7554/eLife.88334.3

## TaG-EM Barcode Analysis Script
TaG-EM analysis scripts and files

### Prerequisites
Python 3

BioPython

cutadapt

### Usage
usage: TaG-EM_barcode_analysis.py [-h] [-i] [-r] [-l] [-o] [-m]

TaG-EM barcode counting script (v1.0) by Daryl Gohl This program takes in a
FASTQ file and barcode reference FASTA file and outputs a plot and a .txt file
of barcode counts.

options:

  -h, --help            show this help message and exit

  -i , --input_folder   Input folder containing raw FASTQ files [required].
  
  -r , --reference_file 
                        Input reference FASTA file [required].
  
  -l , --length         Barcode length (default: 14).
  
  -o , --output_dir     Output directory for barcode count file (default: same
                        folder as input file)
  
  -m , --mismatches_allowed 
                        Number of mismatches to barcode reference sequences
                        allowed (default 2)
### Usage example
TaG-EM_barcode_analysis.py -i <PathToFile/InputFileName> -r <PathToFile/TagEMReferenceFileName>


## TaG-EM Barcode/Cell Barcode Correlation Script
TaG-EM barcode and 10x cell barcode correlation script

### Prerequisites
Python 3

BioPython

### Usage
usage: TaG-EM_barcode_Cell_barcode_correlation.py [-h] [-r1] [-r2] [-x] [-w]
                                                  [-o] [-s] [-m]

TaG-EM barcode and 10x cell barcode correlation script (v1.0) by Daryl Gohl
This program takes in read 1 and read 2 FASTQ files from an TaG-EM barcode
enriched 10x Genomics library, as well as a filtered metadata .csv file and a
barcode reference FASTA file (whitelist) and outputs an updated .csv metadata
file that associates TaG-EM barcode and UMI counts with cell barcodes.

options:

  -h, --help            show this help message and exit
  
  -r1 , --r1_input_file 
                        Input read 1 FASTQ or FASTQ.gz files from enriched
                        TaG-EM library. For reading in 10x cell barcode and
                        UMI [required].
                        
  -r2 , --r2_input_file 
                        Input read 1 FASTQ or FASTQ.gz files from enriched
                        TaG-EM library. For reading in TaG-EM barcode
                        [required].
                        
  -x , --metadata_file 
                        Filtered scRNA-Seq metadata .csv file [required].
                        
  -w , --TagEM_barcode_file 
                        TaG-EM barcode reference FASTA file
                        (whitelist)[required].
                        
  -o , --output_dir     Output directory for barcode count file (default: Home
                        directory)
                        
  -s , --subsample      Number of subsampled reads to analyze (default
                        1000000)
                        
  -m , --mismatches_allowed 
                        Number of mismatches to barcode reference sequences
                        allowed (default 1)

### Usage example
TaG-EM_barcode_Cell_barcode_correlation.py -r1 <PathToFile/Read1FastqFile> -r2 <PathToFile/Read2FastqFile> -x <PathToFile/metadatacsvFile> -o <PathToFile/OutputDirectory> -w <PathToFile/TagEMReferenceFileName> -s <DesiredSubsamplingDepth>
<PathToFile/ReferenceFileName>


## TaG-EM Barcode Reference Files
### Original TaG-EM barcode reference file (initial 20 lines)
TaG-EM_barcodes_v1.fasta

### Extended TaG-EM barcode reference file (176 lines)
TaG-EM_barcodes_v2_extended.fasta
