#!/usr/bin/env python3
import os
import sys
import argparse
import glob
import re

Usage = """
pipeline to clean raw reads from fastq files using HTSstream.
Note: HTSstream has been found to struggle with completely removing adapter sequencing from raw reads of historical sequence data.
SampID is based on a regex search for illumina sequence info from novaseq, newer illumina may require editing to this
part of the script.
Last updated 14 September 2023 Phred M. Benham
"""

##########################################################################################
def HTS_pipe(s, R1, R2,adapters,outdir):

	#include call to HTSstream in a definition, then call definition each time. 
	"""
	path_to_log_file = 
	SampID = 
	R1 = R1 forward sequence
	R2 = R2 reverse sequence
	"""
	
	#output directory
	#within output directory if log_out does not exist mkdir log_out
	log_out = outdir + "stats_log_out"
	if not os.path.exists(log_out):
		os.mkdir(log_out)
	#within output directory if clean_out does not exit mkdir clean_out
	clean_out = outdir + "CleanedReads_out"
	if not os.path.exists(clean_out):
		os.mkdir(clean_out)

	cmd="hts_Stats -L " + log_out + "/" + s + "_stats.log -1 " + R1 + " -2 " + R2 + " -N 'inital raw data stats' | "
	cmd+="hts_SuperDeduper -A " + log_out + "/" + s + "_stats.log -N 'remove duplicates' | "
	cmd+="hts_Overlapper -p 4 -A " + log_out + "/" + s + "_stats.log -N 'trim adaptors' | "
	cmd+="hts_SeqScreener -s " + adapters + " -A " + log_out + "/" + s + "_stats.log -N 'trim any remaining adaptors' | "
	cmd+="hts_QWindowTrim -A " + log_out + "/" + s + "_stats.log -N 'trim poor quality ends' | "
	cmd+="hts_LengthFilter -m 35 -A " + log_out + "/" + s + "_stats.log -N 'remove reads shorter than 35bp' | "
	cmd+="hts_NTrimmer -e -A "  + log_out + "/" + s + "_stats.log -N 'remove any remaining sequences with Ns' | "
	cmd+="hts_Stats -A " + log_out + "/" + s + "_stats.log -F -u -f " + clean_out + "/" + s + " -N 'final stats'" 
	
	return cmd

##########################################################################################

parser = argparse.ArgumentParser(description = Usage)
parser.add_argument("-i","--input", required=True, help="input directory with raw reads")
parser.add_argument("-o","--output", required=True, help="output final cleamed fq files")
parser.add_argument("-a","--adapters", required=True, help="file with adapter sequence to be removed")
args = parser.parse_args()

rawreads = args.input
outdir = args.output
adapters = args.adapters

path = rawreads + "*_R1*.gz"
#print(path)
RawFiles = glob.glob(path)

for r1 in RawFiles:
	#print(r1)
	r2 = r1.replace("R1","R2")
	sampFull = r1.split('/')[-1]
	sampID = re.sub(r'\_S\d+\_L\d+\_.+', '',sampFull)
	print(sampID)
	os_cmd = HTS_pipe(sampID, r1, r2, adapters, outdir)
	os.system(os_cmd)
