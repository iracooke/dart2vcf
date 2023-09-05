#!/usr/bin/env python3

import argparse
import sys
import math
import itertools
import re
import logging
import subprocess
from operator import itemgetter

log = logging.getLogger()
logging.basicConfig(stream=sys.stderr,level=logging.DEBUG)

parser = argparse.ArgumentParser(description='Convert DArT 2 line cvs to vcf')

parser.add_argument('incsv',metavar='FILE',nargs='?',type=argparse.FileType('r'),help='DaRT csv file',default=sys.stdin)

parser.add_argument('-g', '--genome',type=str,help='Genome fasta. Used to find genomic coordinates',required=False)


args  = parser.parse_args()


# First read the header
n_skip=0
samples_start_col=0
info_names = None

for line in args.incsv:
	line_values = line.strip().split(",")

	if line_values[0]=="*":
		n_skip+=1
		if samples_start_col == 0:
			for c in line_values:
				if c=="*":
					samples_start_col+=1
				else:
					break;
	else:
		info_names=line_values
		break;

log.info("Skipped "+str(n_skip)+" header lines. Samples start at column "+str(samples_start_col))


# Look for column indexes of important columns
allele_sequence_c=None
snp_id_c=None

snp_count_c=None
ref_count_c=None

ci=0
for colname in info_names:

	if colname=="AlleleSequence":
		allele_sequence_c=ci

	if colname=="MarkerName" or colname=="AlleleID":
		snp_id_c=ci

	if colname=="AvgCountRefRaw" or colname=="AvgCountRef":
		ref_count_c=ci

	if colname=="AvgCountSnpRaw" or colname=="AvgCountSnp":
		snp_count_c=ci

	ci+=1


# If we specified the `-g` option, open a file for writing the sequences
seqs_file = None

if args.genome is not None:
	seqs_filename="seqs.fasta"
	log.info("Specified a genome. Will write sequences to "+seqs_filename)
	seqs_file = open(seqs_filename,'w')


sample_names = []

for s in info_names[samples_start_col:]:
	suff=""
	ndups=sample_names.count(s)
	if ndups>0:
		suff="_"+str(ndups)
	sample_names.append(s+suff)



n_samples=len(info_names[samples_start_col:])

vcf_records={}

allele_re = re.compile("([ATGC])>([ATGC])")



for line1,line2 in itertools.zip_longest(*[args.incsv]*2):
	line1_values = line1.strip().split(",")
	line2_values = line2.strip().split(",")

	CHROM="chr0"
	POS="0"
	ID=line1_values[snp_id_c]


	m = allele_re.search(ID)
	allele_maj = m.group(1)
	allele_min = m.group(2)

	REF=allele_maj
	ALT=allele_min

	QUAL="."
	FILTER="."

	ref_count = int(float(line2_values[ref_count_c]))
	snp_count = int(float(line2_values[snp_count_c]))

	INFO="DP="+str((ref_count+snp_count)*n_samples)

	FORMAT="GT:AD"

	geno_1=line1_values[samples_start_col:]
	geno_2=line2_values[samples_start_col:]



	genotypes = [ a1.replace('-','.')+"/"+a2.replace('-','.')+":"+str(ref_count)+","+str(snp_count) for a1,a2 in zip(geno_1,geno_2) ]

	vcf_row = [CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT]
	vcf_row.extend(genotypes)

	vcf_records[ID] = vcf_row

	if seqs_file:
		seqs_file.write(">"+ID+"\n"+line1_values[allele_sequence_c]+"\n")


# Map sequences
log.info("Mapping sequences with bwa")

seqs_bam_filename="seqs.bam"

#bwa_args=["bwa","mem","-o",seqs_bam_filename,args.genome,seqs_filename]
bwa_args=["bwa","mem",args.genome,seqs_filename]

bwa_result = subprocess.Popen(bwa_args,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

mapping_records={}

for line in iter(bwa_result.stdout.readline, b""):
	line_parts = line.decode('ascii').strip().split('\t')
	mapping_records[line_parts[0]] = mapping_records.get(line_parts[0],[])
	mapping_records[line_parts[0]].append(line_parts)

for line in iter(bwa_result.stderr.readline, b""):
	log.info(line.decode('ascii').strip())


n_dups=0
n_unmapped=0
n_unique=0
n_indels=0

perfect_cigar_re = re.compile("^([0-9]+)M$")
snp_pos_re = re.compile("([0-9]+):[AGTC]>[AGTC]")

# DArT IDs look like this 26525701|F|0-13:T>A-13
# Note that the numerical range 0-13 indicates that there are 13 bases before the SNP
# In other words, the SNP position is at 14 in a 1-based system or 13 in 0-based
# VCF uses 1-based coordinates

retained_vcfrows=[]

for vcfid,vcfrow in vcf_records.items():
	mappings = mapping_records.get(vcfid,[])
	if len(mappings)==0:
		n_unmapped+=1
	elif len(mappings)>1:
		n_dups+=1
	elif perfect_cigar_re.search(mappings[0][5]) is None:
		n_indels+=1
	else:

		flag=int(mappings[0][1])
		# 4 : Read unmapped
		# 16 : Read reverse strand
		# 256 : Not primary alignment
		# 
		if (flag & 4) or (flag & 256):
			n_unmapped+=1
		else:
			vcfrow[0] = mappings[0][2]
			base_pos = int(mappings[0][3]) # SAM pos so 
			snp_pos_m = snp_pos_re.search(vcfid)
			snp_pos = None
			if snp_pos_m is None:
				log.error("Unable to parse SNP position from ID "+vcfid)
			else:
				snp_pos=int(snp_pos_m.group(1)) # See note on DArT IDs above

			seqlen = len(mappings[0][9])
			if (flag & 16):
				vcfrow[1]=str(base_pos + seqlen - snp_pos -1)
			else:
				vcfrow[1]=str(base_pos + snp_pos)

			retained_vcfrows.append(vcfrow)


log.info("Retained "+
	str(len(retained_vcfrows))+
	" of "+
	str(len(vcf_records))+
	" loci. "+
	str(n_dups)+" with duplicate mappings. "+
	str(n_unmapped)+" unmapped. "+
	str(n_indels)+" with complex CIGAR strings")



# Writing the vcf
vcfout=sys.stdout


# Construct the vcf header

vcfout.write('##fileformat=VCFv4.3'+'\n')
vcfout.write('##source=dart2vcf'+'\n')
vcfout.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate Total Depth">'+'\n')
vcfout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'+'\n')
vcfout.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed. Only locus average is provided">'+'\n')

vcf_samples_row = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
vcf_samples_row.extend(sample_names)

vcfout.write('\t'.join(vcf_samples_row)+'\n')

# import pdb; pdb.set_trace()
# print("hi")

s1=sorted(retained_vcfrows,key=lambda row: int(row[1]))
s2=sorted(s1, key=itemgetter(0))


for vcfrow in s2:
	vcfout.write('\t'.join(vcfrow)+'\n')



