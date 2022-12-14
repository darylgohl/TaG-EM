# TaG-EM Barcode Analysis
TaG-EM analysis scripts and files

## Prerequisites
Python 3

BioPython

cutadapt

## Usage
usage: TaG-EM_barcode_analysis.py [-h] [-i] [-r] [-l] [-o] [-m]

TaG-EM barcode counting script (v1.0) by Daryl Gohl This program takes in a
FASTQ file and barcode reference FASTA file and outputs a plot and a .txt file
of barcode counts.

optional arguments:

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
## Usage example
TaG-EM_barcode_analysis.py -i <PathToFile/InputFileName> -r <PathToFile/ReferenceFileName>
