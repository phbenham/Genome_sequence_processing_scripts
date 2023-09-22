import sys,os,subprocess
import numpy as np
from pyfaidx import Fasta
import multiprocessing as mp
import argparse
import glob
from Bio import SeqIO

##########################################################################################
#get length of full genome and length divided by number of threads to parallelize
def get_fasta_regions(filename, threads):	
	refseq = Fasta(filename)

	key = list(refseq.keys())
	regions = []
	small_scaff = []
	for i in range(len(key)):
		if i < threads-1:
			regions.append(key[i])
	
		else:
			small_scaff.append(key[i])

	regions.append(small_scaff)
	
	return regions

def worker(region,out,bam_files):
    """Run bcftools for a single region."""

    cmd = 'bcftools mpileup -Ob -d 1000 -r {reg} -f {r} -b {b} -a SP,DP -o {o} --threads {t}'.format(r=ref, reg=region, b=bam_files, o=out, t=threads)
    #print (cmd)
    subprocess.check_output(cmd, shell=True)
    cmd = 'bcftools index {o}'.format(o=out)
    subprocess.check_output(cmd, shell=True)
    return
	
def mpileup_parallel(bam_files, ref, outpath, threads, callback=None):
	"""Run mpileup in parallel over multiple regions, then concat vcf files."""

	bam_files = bam_files #' '.join(bam_files) 
	#, replace with bam file list as input
	#bam_file_list = bam_file_list
	rawbcf = os.path.join(outpath,'raw.bcf')
	tmpdir = outpath + 'tmp_mpileup'
	
	pool = mp.Pool(threads)    
	outfiles = []    
	
	
	#get chromosome/scaffold names from get_fasta_regions definition above
	fasta_blocks = get_fasta_regions(ref,threads)
	block_names = []
	for n in range(threads):
		block = "block" + str(n)
		block_names.append(block)

	block_dict = dict(zip(block_names,fasta_blocks))
	        
	for key in block_dict:        
		print (key)
		if len(block_dict[key])<12:
			region = block_dict[key] 
		else:
			region = ','.join(block_dict[key])
			out = '{o}/{s}.bcf'.format(o=tmpdir,s=key)
			f = pool.apply_async(worker, [region,out,bam_files])
			outfiles.append(out)

	pool.close()
	pool.join()

	#concat files
	cmd = 'bcftools concat {i} -O b -o {o}'.format(i=' '.join(outfiles),o=rawbcf)
	subprocess.check_output(cmd, shell=True)
	#remove temp files
	for f in outfiles:
		os.remove(f)
	return rawbcf
##########################################################################################

Usage = "multi-threading bcftools"

parser = argparse.ArgumentParser(description = Usage)
parser.add_argument("-b","--bam", required=True, help="input bam file directory")
parser.add_argument("-f","--fasta", required=True, help="input reference fasta file")
parser.add_argument("-o","--out_path", required=True, help="pathway to output dir")
parser.add_argument("-t","--threads", required=True, help="number of threads")
args = parser.parse_args()



#input options to be user specified
bam_files = args.bam
ref = args.fasta
outpath = args.out_path
threads = int(args.threads)


#use glob to get all bam files from directory
path_name = bam_files + "*.bam"
bamfiles = glob.glob(path_name)

#mkdir tmp_mpileup
if not os.path.exists(outpath + "tmp_mpileup"):
	os.mkdir(outpath + "tmp_mpileup")
	
#mkdir fasta_out
out_fasta = outpath + "fasta_out"
if not os.path.exists(out_fasta):
	os.mkdir(out_fasta)
	
rawbcf = mpileup_parallel(bamfiles, ref, outpath, threads)
	
call_cmd = "bcftools call -m -Ou " + outpath + "raw.bcf | bcftools norm -f " + ref + " -Ou | bcftools filter --IndelGap 5 -i 'QUAL>20 && INFO/DP>10' -Oz -o " + outpath + indv + ".calls.norm.filter.vcf.gz"
os.system(call_cmd)
