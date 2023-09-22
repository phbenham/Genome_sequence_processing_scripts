# Genome_sequence_processing_scripts (last update: 22 September 2023)

This directory contains a series of python scripts for processing raw sequencing
data from WGS or exome capture type projects.

There are four main steps to this pipeline that will clean reads using HTStream, map cleaned reads to
a reference with BWA-MEM, provide basic summary stats on the sequence data, and perform genotype calls with bcftools 

**(1) CleanReads_HTSstream**

	basic usage:
	python3 CleanReads_HTSstream.py -i directory/rawreads/ -o cleanedreads/out/ -a adapters.txt	
		
	Required flags:	
	-i give path to input directory
	-o path to output directory where cleanedreads and json formatted log files will be stored
	-a file with sequence adapters to be trimmed. file adapters.txt has Illumina adapter sequences.	
		
**(2) MapToReference**

	basic usage:
	python3 MapToReference.py -i cleanedreads/out/ -o mapped_reads_out/ -f refgenome.fna -t 12
		
	Required flags:
	-i path to directory with the cleaned fastq files (potentially generated using script above)
	-o path to output directory for mapped reads
	-f reference genome to map reads to
	-t number of threads 
		
**(3) sequence_evaluation_script**

	basic usage:
	python3 sequence_evaluation_script.py -j cleanedreads/out/stats_log_files/ -i mapped_reads_out/ -s sample-data.csv -o evaluation_stats_out/ -n MyEvaluationMetrics
		
	Required flags:
	-j path to directory containing log files in JSON format, output from HTSstream
	-i path to directory containing mapped bam files
	-s File with sample IDs and other meta-data
	-o path to directory for output files
	-n prefix to add to output files
		
	Additional options:
	-b optional bed file for estimating region specific coverage
	-e 1 yes or 0 no, if 1 will calculate specificity and sensitivity for capture datasets [default 0]

**(4) BCFtools_parallel_mpileup**

	basic usage:
	python3 BCFtools_parallel_mpileup.py -b mapped_reads_out/ -f refgenome.fna -o outputvcf/ -t 12

	Required flags:
	-b input bam file directory
	-f input reference fasta file
	-o pathway to output dir
	-t number of threads

**Software dependencies for this pipeline**

 [HTStream](https://s4hts.github.io/HTStream/)
 
 [BWA-MEM](https://bio-bwa.sourceforge.net)
 
 [Sambamba](https://lomereiter.github.io/sambamba/)
 
 [qualimap](http://qualimap.conesalab.org)
 
 [BCFtools](https://samtools.github.io/bcftools/bcftools.html)

python dependencies:

	os
	sys
	argparse
	glob
	re
	json
	numpy
	pandas
	matplotlib
	seaborn
	multiprocessing
	pyfaidx
	Bio

If installing on a Linux machine the CCGPpipeline.yml can be used to install all
dependencies with the following command (assuming you have miniconda or anaconda installed):

conda env create -f CCGPpipeline.yml
