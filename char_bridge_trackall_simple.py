#!/usr/bin/python
USAGE="""Nikki Teran 2016-03-21 from program started 2016-1-16
python char_bridge.py --R1 R1.fastq.gz --R2 R1.fastq.gz
python char_bridge.py --ASSEMBLED S0.assembled.fastq.gz

Separates single end reads into RNA and DNA based on an orientation specific bridge.

Input files must be gzipped
"""
import os
import re
import sys
import bz2
import time
import gzip
import string
import argparse
import charbridgetools as cbt
from Bio import SeqIO
from Bio import AlignIO
from itertools import izip
from optparse import OptionParser

#first things first
bridge1= 'AAACCGGCGTCCAAG'  # 'AACATAAACCGGCGTCCAAGGATCTTTAATTAAGTCGCAG'  # before ligation with last ATC 'AACATAAACCGGCGTCCAAGGATCTTTAATTAAGTCGCAGATC' # after ligation = 'AACATAAACCGGCGTCCAAGGATC'  # default = 'ACCGGCGTCCAAG' #fwd, so RNA, bridge, DNA
#bridge2='CTTGGACGCCGGT' #rev, so DNA, bridge, RNA
b1_len=len(bridge1)
#b2_len=len(bridge2)

#all the arguments we want to parse.
# python char_bridge_trackall.py --FASTQGZ meaningless.fasta.gz --NAME test. --minRNA 2 --minDNA 2
parser = argparse.ArgumentParser(description = USAGE)
parser.add_argument('--FASTQGZ', dest = 'FASTQGZ', default = "R1.fastq.gz", help = 'REQUIRED FOR PAIRED END READ FUNCTIONALITY: Read 1 input file in gzipped fastq format or forward unassembled read from PEAR')

parser.add_argument('--DNA', dest = 'DNA', default = "dna.fastq.gz", help = 'Desired output for the DNA sequences')
parser.add_argument('--RNA', dest = 'RNA', default = "rna.fastq.gz", help = 'Desired output for the RNA sequences')

parser.add_argument('--NB', dest = 'NB', default = ".nobridge.fastq.gz", help = 'Desired output for Read 1 lines that do not contain the bridge')

parser.add_argument('--DB', dest = 'DB', default = "dupbridge.fastq.gz", help = 'Desired output for Read 1 lines that contain duplicate bridges')

parser.add_argument('--TS', dest = 'TS', default = "tooshort.fastq.gz", help = 'Desired output for RNA when either the RNA or DNA is too short')

parser.add_argument('--POS', dest = 'POS', default = "bridgeposition.fastq.gz", help = 'Desired output for list of bridge positions, if bridge is found')

##### VCM added this code below #####
parser.add_argument('--TL', dest = 'TL', default = "RNAtoolong.fastq.gz", help = 'Desired output for RNA when either the RNA or DNA is too short')
##### VCM added this code above #####

parser.add_argument('--minRNA', dest = 'minRNA', default = 18, help = 'The minimum length of RNA, shorter reads will be discarded. Default: 18')
parser.add_argument('--minDNA', dest = 'minDNA', default = 18, help = 'The minimum length of RNA, shorter reads will be discarded. Default: 18')
##### VCM added this code below #####
parser.add_argument('--maxRNA', dest = 'maxRNA', default = 1000, help = 'The minimum length of RNA, shorter reads will be discarded. Default: 18')
##### VCM added this code above #####

parser.add_argument('--NAME', dest = 'NAME', default = "", help = 'Sample name to add to start of all output files.')

args = parser.parse_args()
name=args.NAME

R1_fh = args.FASTQGZ
DNA_fh = name+args.DNA
RNA_fh = name+args.RNA
POS_fh = name+args.POS
##### VCM added this code below #####
AANNN_fh = name+'AANNN.gz'
##### VCM added this code above #####


minDNA = args.minDNA
minRNA = args.minRNA
##### VCM added this code below #####
maxRNA = args.maxRNA
##### VCM added this code above #####


#open files
rna_out = gzip.open(RNA_fh, 'w')
dna_out = gzip.open(DNA_fh, 'w')
pos_out = gzip.open(POS_fh, 'w')
##### VCM added this code below #####
aannn_out = gzip.open(AANNN_fh, 'w')
##### VCM added this code above #####

pos_temp, aannn_temp, rna_temp, dna_temp = [], [], [], []

with gzip.open(R1_fh) as read_in:
	linecount=0
	#for line in read_in:
	for count, line in enumerate(read_in, start=0):
		line=line.strip()
		linecount=linecount+1
		if linecount==1:  # if line =1 save the name
			name=line
		elif linecount==2:  # if line =2, hunt for the bridge
			ntseq=line
			rnaseq,dnaseq,orientation,position,bridgenum,NNN,nnngate,dnagate=cbt.bridgehunter(bridge1,line)  # VCM added
		elif linecount==4:  # if line =3 plus sign # if line == 4 split the quality score
			linecount=0
			if orientation=='R':  #remember to flip the line if we found the string backwards
				line=line[::-1]
			rnaqual=line[0:len(rnaseq)]
			dnaqual=line[-len(dnaseq):]
			#print out the position information for everything. Should be as 1/4 as many lines in pos as in input file
			pos_temp.append(str(position)+'\t'+str(len(line)-position-len(bridge1))+'\n')
			if (bridgenum == 1) and (nnngate == 1) and (dnagate == 1) and (len(dnaseq) >= int(minDNA)) and (len(rnaseq) >= int(minRNA)) and (len(rnaseq) <= int(maxRNA)):
				# outname = name  #VCM added NNN to seq name for PCR duplicate removal
				seqname = ':'.join([name.split()[0], name.split()[1].split(':')[0], NNN]) + ' ' + name.split()[1]
				rna_temp.append('\n'.join([seqname, rnaseq, '+', rnaqual]))
				dna_temp.append('\n'.join([seqname, dnaseq, '+', dnaqual]))
				aannn_temp.append(NNN)
			else:
				pass
		if count % 12000 == 0:  # output to file every 3000 seqs
			if len(pos_temp) > 0:
				pos_out.write('\n'.join(pos_temp) + '\n')
			if len(aannn_temp) > 0:  # VCM added #
				aannn_out.write('\n'.join(aannn_temp) + '\n')  # VCM added #
			if len(rna_temp) > 0:
				rna_out.write('\n'.join(rna_temp) + '\n')
			if len(dna_temp) > 0:
				dna_out.write('\n'.join(dna_temp) + '\n')
			pos_temp, aannn_temp, rna_temp, dna_temp = [], [], [], []

	if len(pos_temp) > 0:
		pos_out.write('\n'.join(pos_temp) + '\n')
	if len(aannn_temp) > 0:  # VCM added #
		aannn_out.write('\n'.join(aannn_temp) + '\n')  # VCM added #
	if len(rna_temp) > 0:
		rna_out.write('\n'.join(rna_temp) + '\n')
	if len(dna_temp) > 0:
		dna_out.write('\n'.join(dna_temp) + '\n')

##### VCM added this code below #####
#close files
rna_out.close()
dna_out.close()
aannn_out.close()
pos_out.close()
##### VCM added this code above #####
