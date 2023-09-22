#!/usr/bin/env python3
import os
import sys
import argparse
import glob
import re

Usage = """
pipeline to align fastq files to reference genome followed by processing of bam files.
programs used: bwa, sambamba, picard tools, gatk3
Note: input fastq names can be either compressed or uncompressed, but must all have the 
syntax *_R1*, *_R2*, and *_SE*.
version 2.0, updated by Phred M. Benham on 9 May 2023
"""

parser = argparse.ArgumentParser(description = Usage)
parser.add_argument("-i","--input", required=True, help="input cleaned reads")
parser.add_argument("-o","--output", required=True, help="output final mapped bam files")
parser.add_argument("-f","--reference", required=True, help="reference genome")
parser.add_argument("-t","--threads", required=True, help="number of threads")
args = parser.parse_args()

#input options to be user specified on command line
CleanedReads=args.input
MappedData=args.output
RefGenome=args.reference
NumThreads=int(args.threads)


path_name = CleanedReads + "*_R1*"
files = glob.glob(path_name)

#if statement to search for indexed bwa reference
#if it does not exist, move forward, else make index.
bwa_ref_index = RefGenome + ".bwt"
if not os.path.exists(bwa_ref_index):
	os.system("bwa index %s" % RefGenome)

#for loop over each sample 
for file in files:
	Fastq1 = file
	results = re.findall('.+\/(.+)\_R1.+',file)
	sample_name=results[0]
	print(sample_name)
	Fastq2 = re.sub(r'\_R1', '_R2',file)
	FastqU = re.sub(r'\_R1', '_SE',file)
	OutSamPaired =  MappedData + sample_name + ".outPairedSam1"
	OutSamSolo = MappedData + sample_name + ".outSoloSam1"
	OutbamPaired =  MappedData + sample_name + ".outPaired.bam"
	OutbamSolo = MappedData + sample_name + ".outSolo.bam"
	SortedPairedBam = MappedData + sample_name + ".outSortedPaired.bam"
	SortedSoloBam = MappedData + sample_name + ".outSortedSolo.bam"
	MergedSortedBam = MappedData + sample_name + ".outSortedMerged.bam"	
	sorted_in_target_bams = MappedData + sample_name + "_sorted";
	sorted_in_target_bams2 = sorted_in_target_bams + '.bam';

	os.system("bwa mem -t %d %s %s %s > %s" % (NumThreads, RefGenome, Fastq1, Fastq2, OutSamPaired))
	os.system("bwa mem -t %d %s %s > %s" % (NumThreads, RefGenome, FastqU, OutSamSolo))

	#step 3 post-processing of sam/bam files
	
	#remove non-uniquely mapping reads and convert to bam
	os.system("sambamba view -t %d -h -f bam -S -F 'not ([XA] != null or [SA] != null)' %s -o %s" % (NumThreads,OutSamPaired,OutbamPaired))
	os.system("sambamba view -t %d -h -f bam -S -F 'not ([XA] != null or [SA] != null)' %s -o %s" % (NumThreads,OutSamSolo,OutbamSolo))
	
	os.unlink(OutSamPaired)
	os.unlink(OutSamSolo)
	
	#sort bamfiles
	os.system("sambamba sort -t %d -o %s %s" % (NumThreads, SortedPairedBam, OutbamPaired))
	os.system("sambamba sort -t %d -o %s %s" % (NumThreads, SortedSoloBam, OutbamSolo))
	os.unlink(OutbamPaired)
	os.unlink(OutbamSolo)
	#merge bamfiles
	os.system("sambamba merge -t %d %s %s %s" % (NumThreads, MergedSortedBam, SortedPairedBam, SortedSoloBam))
	os.unlink(SortedPairedBam)
	os.unlink(SortedSoloBam)
	#index bamfiles
	os.system("sambamba index -t %d %s" % (NumThreads, MergedSortedBam))

