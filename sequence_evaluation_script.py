#!/usr/bin/env python3
import os
import sys
import argparse
import glob
import re
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('notebook')

Usage = """
Script to calculate a variety of QC metrics for assessing overall quality of data.
Supply input file with sample ID and other meta-data to add on columns. 
Uses log output files from HTSstream to calculate total amount of raw data, cleaned data, etc. outputs figures to aid quality assessment.
version 1.0, created by Phred M. Benham on 23 February 2023
"""

#input flags:
parser = argparse.ArgumentParser(description = Usage)
parser.add_argument("-j","--Json_dir", required=True, help="path to directory containing log files in JSON format, output from HTSstream")
parser.add_argument("-i","--Bam_files_dir", required=True, help="path to directory containing bam files")
parser.add_argument("-s", "--Sample_data", required=True, help="File with sample IDs and other meta-data")
parser.add_argument("-b", "--bed_file", required=False, help="optional bed file for estimating region specific coverage")
parser.add_argument("-e","--estimate_sensitivity_specificity", nargs='?', const=1, type=int, default=0, help="1 yes or 0 no, if 1 will calculate specificity and sensitivity for capture datasets [default 0]")
parser.add_argument("-o","--output_dir", required=True, help="path to directory for output files")
parser.add_argument("-n","--Outfile_prefix", required=True, help="prefix to add to output files")

args = parser.parse_args()
rawlog = args.Json_dir
bams = args.Bam_files_dir
samples = args.Sample_data
bedfile = args.bed_file
SensSpec = args.estimate_sensitivity_specificity
OutDir = args.output_dir
OutName = args.Outfile_prefix

#mk output directory if it does not exist.
if not os.path.exists(OutDir):
	os.mkdir(OutDir)

#import samples file
SampleFile = pd.read_csv(samples, index_col=0)
print(SampleFile.head())

path = rawlog+"*.log"
RawFiles = glob.glob(path)

ColNames1 = ["SampID","StartReads","DeDuplicate", "OverLapper", "AdapterRemoval", "QualityTrim", "LengthFilter","NTrimmer","EndBP", "EndReads", "StartBP"]
#delete ColNames = ["SampID","DeDuplicate", "OverLapper", "AdapterRemoval", "QualityTrim", "LengthFilter","NTrimmer"]

for newCol in ColNames1 :
	SampleFile[newCol] = ''

c=0
for log in RawFiles:
	indv=log.split('/')[-1].split('_')[0]
	print(indv)
	SampleFile.loc[indv,"SampID"] = indv
	
	#parse json log file from HTS stream
	logfile_df = pd.read_json(log)
	logfile_df = logfile_df.drop(labels=[0,7])
	fragments = logfile_df['Fragment']
	Start_reads = fragments[1]['in']
	Start_bp = fragments[1]['basepairs_in']
	SampleFile.loc[indv,"StartReads"] = Start_reads
	SampleFile.loc[indv,"StartBP"] = Start_bp
	
	n=1
	for row in fragments:
		SampleFile.loc[indv,"EndReads"] = row['out']
		SampleFile.loc[indv,"EndBP"] = row['basepairs_out']
		SampleFile.loc[indv,ColNames1[n+1]] = row['basepairs_out']
		n+=1	

SampleFile["Percent_reads_retained"] = (SampleFile["EndReads"]/SampleFile["StartReads"])*100
SampleFile["Percent_bp_retained"] = (SampleFile["EndBP"]/SampleFile["StartBP"])*100

#perform multi-bamqc using qualimap

# first need to create file for -d flag in command line
BamPath = bams+"*.bam"
RawBam = glob.glob(BamPath)

#mkdir for qualimap_results
qualmap_dir = bams + "qualimap_out"
if not os.path.exists(qualmap_dir):
        os.mkdir(qualmap_dir)

BamList = bams + OutName + "_bamlist.txt"      
OutFile = open(BamList, 'w') 
for bam in RawBam:
	Outfilename = bam.split('/')[-1].split('.')[0]
	print(Outfilename, bam, sep='\t', file=OutFile)

OutFile.close()

print(Outfilename)
cmd = "qualimap multi-bamqc -d " + BamList + " -gff " + bedfile + " -r  -outdir " + qualmap_dir + " --java-mem-size=4G"
print(cmd)
#os.system(cmd)


NewColNames = ["MeanCov","std_cov", "GC_content", "meanMQ", "PCTreads_mapped_OnTarget", "PercentTarget_1xCov","2x_cov", "5x_cov", "10x_cov"]
for newCol1 in NewColNames:
	SampleFile[newCol1] = ''

#append results to data frame
for bam in RawBam:
	SampID = bam.split('/')[-1].split('.')[0]
	middle= bam.split('/')[-1].split('.')[1]
	#replace .bam with _stats
	resultsFile = bams + SampID + "." + middle + "_stats/genome_results.txt"
	qc_results = open(resultsFile, 'r')
	for line in qc_results:
		line = line.strip('\n')
		
		if "number of mapped reads =" in line:
			PCT_OnTargetReads = float(line.split('(')[1].split(')')[0].split('%')[0])
			SampleFile.loc[SampID,"PCTreads_mapped_OnTarget"] = PCT_OnTargetReads
		
		elif "mean mapping quality =" in line:	
			meanMQ = float(line.split(' ')[-1])
			SampleFile.loc[SampID,"meanMQ"] = meanMQ
		
		elif "GC percentage = " in line:
			GC_pct = float(line.split(' ')[-1].split('%')[0])
			SampleFile.loc[SampID,"GC_content"] = GC_pct
			
		elif "mean coverageData = " in line:
			MeanCov = float(line.split(' ')[-1].split('X')[0])
			SampleFile.loc[SampID,"MeanCov"] = MeanCov
			
		elif "std coverageData = " in line:
			StdCov = float(line.split(' ')[-1].split('X')[0])
			SampleFile.loc[SampID,"std_cov"] = StdCov
			
		elif " >= " in line:
			CovLevel = line.split(' ')[-1]
			if CovLevel == "1X":
				PctTargets_1x_cov = re.findall('(\d*\.\d*)%',line)[0]
				SampleFile.loc[SampID,"PercentTarget_1xCov"] = PctTargets_1x_cov
				
			elif CovLevel == "2X":
				PctTargets_2x_cov = re.findall('(\d*\.\d*)%',line)[0]
				SampleFile.loc[SampID,"2x_cov"] = PctTargets_2x_cov
				
			elif CovLevel == "5X":
				PctTargets_5x_cov = re.findall('(\d*\.\d*)%',line)[0]
				SampleFile.loc[SampID,"5x_cov"] = PctTargets_5x_cov
				
			elif CovLevel == "10X":
				PctTargets_10x_cov = re.findall('(\d*\.\d*)%',line)[0]
				SampleFile.loc[SampID,"10x_cov"] = PctTargets_10x_cov			

SampleFileName = OutDir + OutName + "_Out.csv"
SampleFile.to_csv(SampleFileName)


#make figures
#make figure of percent bp retained individual
fig, ax = plt.subplots(figsize=(12, 6))
sns.despine(ax=ax, offset=5)
ax = sns.barplot(x='SampID', y='Percent_bp_retained',hue="Tissue",  data=SampleFile, order=SampleFile.sort_values('Percent_bp_retained').SampID, dodge=False)
ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)	
ax.axhline(SampleFile['Percent_bp_retained'].mean(), color="red", lw=0.5)
ax.set_ylabel("Percent bp retained after cleaning reads")
ax.set(title="Per individuals base pairs(%) retained after cleaning in HTSstream")
fig.tight_layout()
OutFileName = OutDir + OutName + "_percent_retained_bp.pdf"
fig.savefig(OutFileName, dpi=300)

#make figure of total bp retained for each individual
fig, ax = plt.subplots(figsize=(12, 6))
sns.despine(ax=ax, offset=5)
ax = sns.barplot(x='SampID', y="EndBP",hue="Tissue", data=SampleFile, order=SampleFile.sort_values('EndBP').SampID, dodge=False)
ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)	
ax.axhline(SampleFile['EndBP'].mean(), color="red", lw=0.5)
ax.set_ylabel("Total bp retained after cleaning reads (1e10)")
ax.set(title="Total bp retained after cleaning in HTSstream")
fig.tight_layout()
OutFileName2 = OutDir + OutName + "_total_retained_bp.pdf"
fig.savefig(OutFileName2, dpi=300)

#Create boxplot showing retention of bps across different steps
#transpose df then melt
#subset df for melting
#cols to retain ["DeDuplicate", "OverLapper", "AdapterRemoval", "QualityTrim", "LengthFilter", "NTrimmer"]
NewDf = SampleFile[["DeDuplicate", "OverLapper", "AdapterRemoval", "QualityTrim", "LengthFilter", "NTrimmer"]]

#df_per_bp_T = df_per_bp.drop(labels="SampID", axis=1)
NewDf_melt = NewDf.melt()

NewDf_melt['value'] = NewDf_melt['value'].astype(int)

fig, ax = plt.subplots(figsize=(12, 6))
sns.despine(ax=ax, offset=5)
ax = sns.boxplot(x='variable', y='value', data=NewDf_melt, color="gray")
#ax.set_xticklabels(ax.get_xticklabels(), rotation=90)	
#ax.axhline(df_bp['NTrimmer'].mean(), color="red", lw=0.5)
ax.set_ylabel("Bp retained after each cleaning")
ax.set(title="Bp retained after each step")
OutFileName3 = OutDir + OutName + "_step_retention_bp.pdf"
fig.tight_layout()
fig.savefig(OutFileName3, dpi=300)