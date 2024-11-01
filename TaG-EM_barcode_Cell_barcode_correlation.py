#!/usr/bin/env python

"""
Copyright 2024  Daryl Gohl

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from Bio import SeqIO
import csv 
from collections import defaultdict
import gzip
import regex
import os
import argparse
__version__ = "1.0"


import numpy as np
import matplotlib.pyplot as plt

def get_args(x):
    x.add_argument("-r1", "--r1_input_file",
                   type = str,
                   default = None,
                   metavar = '',
                   help="Input read 1 FASTQ or FASTQ.gz files from enriched TaG-EM library. For reading in 10x cell barcode and UMI [required].")
    x.add_argument("-r2", "--r2_input_file",
                   type = str,
                   default = None,
                   metavar = '',
                   help="Input read 1 FASTQ or FASTQ.gz files from enriched TaG-EM library. For reading in TaG-EM barcode [required].")
    x.add_argument("-x", "--metadata_file",
                   type = str,
                   default = None,
                   metavar = '',
                   help="Filtered scRNA-Seq metadata .csv file [required].")
    x.add_argument("-w", "--TagEM_barcode_file",
                   type = str,
                   default = None,
                   metavar = '',
                   help="TaG-EM barcode reference FASTA file (whitelist)[required].")
    x.add_argument("-o", "--output_dir",
                   type = str,
                   default = '',
                   metavar = '',
                   help = "Output directory for barcode count file (default: Current directory)") 
    x.add_argument("-s", "--subsample",
                   type = str,
                   default = '1000000',
                   metavar = '',
                   help = "Number of subsampled reads to analyze (default 1000000)")  
    x.add_argument("-m", "--mismatches_allowed",
                   type = str,
                   default = '1',
                   metavar = '',
                   help = "Number of mismatches to barcode reference sequences allowed (default 1)")  
    args = x.parse_args()
    return args

# Parses command line arguments
argparser = argparse.ArgumentParser(description = "TaG-EM barcode and 10x cell barcode correlation script (v" + __version__ +")\n" + \
                                              "by Daryl Gohl\n" + \
                                              "This program takes in read 1 and read 2 FASTQ files from an TaG-EM barcode enriched 10x Genomics library, as well as a filtered metadata .csv file and a barcode reference FASTA file (whitelist) and outputs an updated .csv metadata file that associates TaG-EM barcode and UMI counts with cell barcodes.",
                                add_help = True, 
                                epilog ='')
args = get_args(argparser)
args = vars(args)

# possible input args
r1_file = args['r1_input_file']
r2_file = args['r2_input_file']
cell_ID_file = args['metadata_file']
barcode_whitelist = args['TagEM_barcode_file']
mismatches = args['mismatches_allowed']
number_of_reads_to_examine = int(args['subsample'])
out_folder = args['output_dir']

if out_folder == '':
  out_dir = os.path.expanduser("~")
else:
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    out_dir = out_folder
    
#Read in list and get cell IDs
columns = defaultdict(list)
with open(cell_ID_file) as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        for (i,v) in enumerate(row):
            columns[i].append(v)
            
#columns[0] = cell IDs
cell_barcode_list = []
for i in columns[0]:
	cell_barcode_list.append(i[0:16])

#Read in sequencing data
#read 1 - get cell barcode and UMI (separately)
R1_cell_barcode_list = []
R1_UMI_list = []
read_count = 0
if r1_file[-6:] == ".fastq":
	for record in SeqIO.parse(handle, "fastq"):
		if read_count < number_of_reads_to_examine:
			xx = str(record.seq)
			R1_cell_barcode_list.append(xx[0:16])
			R1_UMI_list.append(xx[16:26])
			read_count += 1
		else:
			break
elif r1_file[-9:] == ".fastq.gz":
	with gzip.open(r1_file, "rt") as handle:
		for record in SeqIO.parse(handle, "fastq"):
			if read_count < number_of_reads_to_examine:
				xx = str(record.seq)
				R1_cell_barcode_list.append(xx[0:16])
				R1_UMI_list.append(xx[16:26])
				read_count += 1
			else:
				break
		
#read 2 - get TaG-EM barcode
R2_TagEM_barcode_list = []
read_count = 0
if r2_file[-6:] == ".fastq":
    for record in SeqIO.parse(handle, "fastq"):
    	if read_count < number_of_reads_to_examine:
    		read_count += 1
    		xx = str(record.seq)
    		xx_pos = xx.find('CAACAACCGGAAGTGA')
    		if xx_pos != -1:
    			R2_TagEM_barcode_list.append(xx[(xx_pos+16):(xx_pos+30)])
    		else:
    			R2_TagEM_barcode_list.append('Off Target')
    	else:
    		break
elif r1_file[-9:] == ".fastq.gz":
	with gzip.open(r2_file, "rt") as handle:
		for record in SeqIO.parse(handle, "fastq"):
			if read_count < number_of_reads_to_examine:
				read_count += 1
				xx = str(record.seq)
				xx_pos = xx.find('CAACAACCGGAAGTGA')
				if xx_pos != -1:
					R2_TagEM_barcode_list.append(xx[(xx_pos+16):(xx_pos+30)])
				else:
					R2_TagEM_barcode_list.append('Off Target')
			else:
				break

#Get barcode whitelist
TaGEM_seq_list = [] #Barcode sequence
TaGEM_BC_name_list = [] #Barcode name
for record in SeqIO.parse(barcode_whitelist, "fasta"):
    TaGEM_seq_list.append(str(record.seq))
    TaGEM_BC_name_list.append(str(record.id))

#Correlate observed cell barcodes with TaG-EM barcodes
unique_TaGEM_barcode_observed = []
unique_TaGEM_barcode_counts = []
single_unique_TaGEM_barcode = []
cell_barcodes_observed = []
distinct_UMIs = []
for i in cell_barcode_list:
	temp_TaGEM_list = []
	temp_UMI_list = []
	temp_cell_barcode = 'not observed'
	for j, item in enumerate(R1_cell_barcode_list):
		if i == item:
			Whitelist_BC_detected = False
			Whitelist_BC_ID = ''
			for n in TaGEM_seq_list:
				query = r'(?:' + R2_TagEM_barcode_list[j] +'){s<=' + mismatches + '}' #fuzzy matching - allow up to user specified number of mismatches (default = 1)
				test = regex.findall(query, n)
				if test != []:
					Whitelist_BC_ID = n
					Whitelist_BC_detected = True
					temp_TaGEM_list.append(Whitelist_BC_ID)
					temp_UMI_list.append(R1_UMI_list[j])
					temp_cell_barcode = i
	#count unique barcodes
	TaGEM_count = 0
	unique_BC = set(temp_TaGEM_list)
	unique_BC_counts = []
	unique_BC_list = []
	for k in unique_BC:
		if k != 'Off Target':
			TaGEM_count += 1
			temp_count = 0
			for m in temp_TaGEM_list:
				if m == k:
					temp_count += 1
			unique_BC_counts.append(temp_count)
			unique_BC_list.append(k)
	if len(unique_BC_list) == 1:
		single_unique_TaGEM_barcode.append(True)
	else:
		single_unique_TaGEM_barcode.append(False)	
	unique_TaGEM_barcode_observed.append(unique_BC_list)
	unique_TaGEM_barcode_counts.append(unique_BC_counts)		
	distinct_UMIs.append(len(set(temp_UMI_list)))
	cell_barcodes_observed.append(temp_cell_barcode)

#Calculate some summary statistics
most_abundant_TaGEM_barcode = []
most_abundant_TaGEM_barcode_name = []
percent_most_abundant_TaGEM_barcode = []
for i, item in enumerate(unique_TaGEM_barcode_counts):
	if item != []:
		most_abundant_TaGEM_barcode.append(unique_TaGEM_barcode_observed[i][item.index(max(item))])
		percent_most_abundant_TaGEM_barcode.append(max(item)/sum(item)*100)
		for j, item2 in enumerate(TaGEM_seq_list):
			if item2 == unique_TaGEM_barcode_observed[i][item.index(max(item))]:
				most_abundant_TaGEM_barcode_name.append(TaGEM_BC_name_list[j])
	else:
		most_abundant_TaGEM_barcode.append("None")
		percent_most_abundant_TaGEM_barcode.append("NA")
		most_abundant_TaGEM_barcode_name.append("None")
			
#Write out summary_file
outfile_name = 'metadata_with_TagEM_barcodes_' + str(number_of_reads_to_examine) +'.csv'
cell_barcode_corr_file = os.path.join (out_folder, outfile_name)

reader = csv.reader(open(cell_ID_file, 'r'))
writer = csv.writer(open(cell_barcode_corr_file, 'w'))
headers = next(reader)
headers.append("Cell Barcode")
headers.append("TaG-EM Barcodes detected")
headers.append("TaG-EM Barcode counts")
headers.append("UMIs observed")
headers.append("Most abundant TaG-EM barcode name")
headers.append("Most abundant TaG-EM barcode sequence")
headers.append("Percent of most abundant TaG-EM barcode")
writer.writerow(headers)
count = 0
for row in reader:
    row.append(cell_barcodes_observed[count])
    row.append(unique_TaGEM_barcode_observed[count])
    row.append(unique_TaGEM_barcode_counts[count])
    row.append(distinct_UMIs[count])
    row.append(most_abundant_TaGEM_barcode_name[count])
    row.append(most_abundant_TaGEM_barcode[count])
    row.append(percent_most_abundant_TaGEM_barcode[count])
    writer.writerow(row)
    count+=1
    
