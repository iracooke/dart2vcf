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

parser.add_argument('-o', '--vcfout',metavar='FILE',type=argparse.FileType('w'),help='Filename for vcf output',default=sys.stdout)

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


# convert column headers to a dictionary

info_dict = {cn: ci for ci,cn in enumerate(info_names[0:samples_start_col])}


#
# Look for column indexes of essential columns
#
# The ID is usually given by AlleleID but sometimes MarkerName
snp_id_c = info_dict['AlleleID'] if info_dict.get('AlleleID') is not None else info_dict['MarkerName']

# Used for genome mapping
# The trimmed sequence is preferred here since adaptors will lead to clipping or mapping failure
#
allele_sequence_c = info_dict['TrimmedSequence'] if info_dict.get('TrimmedSequence') is not None else info_dict['AlleleSequence']

# Some files have a suffix (raw/norm) appended to this important field
ref_count_c = info_dict['AvgCountRef'] if info_dict.get('AvgCountRef') is not None else info_dict['AvgCountRefRaw']
snp_count_c = info_dict['AvgCountSnp'] if info_dict.get('AvgCountSnp') is not None else info_dict['AvgCountSnpRaw']

# All other INFO columns
# Details on these columns obtained from maptk https://bitbucket.org/rokhsar-lab/gbs-analysis/src/master/bin/maptk
#
#
# The SNP and CloneID fields are not included because they violates formatting requirement for INFO fields. All info in these fields is parsed and included anyway
#
optional_info = {
#	'SNP':['String', 1, '"Contains the base position and base variant details"'],
#	'CloneID':['String', 1, '"Unique identifier of the sequence tag"'],
	'AvgCountRef':['Float', 1, '"The sum of the tag read counts for all samples, divided by the number of samples with non-zero tag read counts, for the Reference allele row"'],
	'AvgCountSnp':['Float', 1, '"The sum of the tag read counts for all samples, divided by the number of samples with non-zero tag read counts, for the SNP allele row"'],
	'AvgPIC':['Float', 1, '"The average of the polymorphism information content (PIC) of the Reference and SNP allele rows"'],
	'CallRate':['Float', 1, '"The proportion of samples for which the genotype call is not absent"'],
	'FreqHets':['Float', 1, '"The proportion of samples which score as heterozygous"'],
	'FreqHomRef':['Float', 1, '"The proportion of samples which score as homozygous for the Reference allele"'],
	'FreqHomSnp':['Float', 1, '"The proportion of samples which score as homozygous for the SNP allele"'],
	'OneRatioRef':['Float', 1, '"The proportion of samples for which the genotype score is 1 in the Reference allele row"'],
	'OneRatioSnp':['Float', 1, '"The proportion of samples for which the genotype score is 1 in the SNP allele row"'],
	'NS':['Integer',1,'"Number of Samples With Data"'],
	'PICSnp':['Float', 1, '"The polymorphism information content (PIC) for the SNP allele row"'],
	'PICRef':['Float', 1, '"The polymorphism information content (PIC) for the Reference allele row"'],
	'RepAvg':['Float', 1, '"The proportion of technical replicate assay pairs for which the marker score is consistent"'],
	'SnpPosition':['Integer', 1, '"The position in the sequence tag at which the defined SNP variant base occurs"'],
	'AlleleSequence':['String', 1, '"The sequence of the Reference allele"'],
	'TrimmedSequence':['String', 1, '"Same as the full sequence, but with removed adapters in short marker tags"']
}

# Column indexes for any of the above columns that exist in the file
optional_info_indexes = { key: info_dict[key] for key in optional_info.keys() if info_dict.get(key)}


# If we specified the `-g` option, open a file for writing the sequences
# This file will be used as input for bwa mem
seqs_file = None

if args.genome is not None:
	seqs_filename="seqs.fasta"
	log.info("Specified a genome. Will write sequences to "+seqs_filename)
	seqs_file = open(seqs_filename,'w')


sample_names = []
sample_tally = []

for s in info_names[samples_start_col:]:
	suff=""
	ndups=sample_tally.count(s)
	if ndups>0:
		suff="_"+str(ndups)
	sample_names.append(s+suff)
	sample_tally.append(s)

#import pdb;pdb.set_trace()


n_samples=len(info_names[samples_start_col:])

vcf_records={}

allele_re = re.compile("([ATGC])>([ATGC])")



for line1,line2 in itertools.zip_longest(*[args.incsv]*2):
	line1_values = line1.strip().split(",")
	line2_values = line2.strip().split(",")

	# Placeholder values for genomic coordinates. If -g is specified these will be updated later
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

	#  (String, no whitespace, semicolons, or equals-signs permitted; commas are permitted only as delimiters for lists of values)
	# 
	INFO = ["DP="+str((ref_count+snp_count)*n_samples)]
	for field,index in optional_info_indexes.items():
		INFO.append(field+"="+str(line2_values[index]))

	INFO = ';'.join(INFO)

	FORMAT="GT:AD"

	geno_1=line1_values[samples_start_col:]
	geno_2=line2_values[samples_start_col:]



	genotypes = [ a1.replace('-','.')+"/"+a2.replace('-','.')+":"+str(ref_count)+","+str(snp_count) for a1,a2 in zip(geno_1,geno_2) ]

	vcf_row = [CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT]
	vcf_row.extend(genotypes)

	vcf_records[ID] = vcf_row

	if seqs_file:
		seqs_file.write(">"+ID+"\n"+line1_values[allele_sequence_c]+"\n")

retained_vcfrows=[]
if args.genome is not None:
	# Map sequences
	log.info("Mapping sequences with bwa")

	seqs_bam_filename="seqs.bam"

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
# Note that the numerical range 0-13 indicates that there are 14 bases before the SNP
# In other words, the SNP position is at 13 in 0-based system
# VCF uses 1-based coordinates



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
else:
	retained_vcfrows=vcf_records.values()

# Writing the vcf
vcfout=args.vcfout


# Construct the vcf header

vcfout.write('##fileformat=VCFv4.3'+'\n')
vcfout.write('##source=dart2vcf'+'\n')
vcfout.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate Total Depth">'+'\n')

for infokey in optional_info_indexes.keys():
	info_type,info_multiplicity,info_desc = optional_info[infokey]
	vcfout.write('##INFO=<ID='+infokey+",Number="+str(info_multiplicity)+",Type="+info_type+",Description="+info_desc+'>'+'\n')


vcfout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'+'\n')
vcfout.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed. Only locus average is provided">'+'\n')

vcf_samples_row = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
vcf_samples_row.extend(sample_names)

vcfout.write('\t'.join(vcf_samples_row)+'\n')

s1=sorted(retained_vcfrows,key=lambda row: int(row[1]))
s2=sorted(s1, key=itemgetter(0))


for vcfrow in s2:
	vcfout.write('\t'.join(vcfrow)+'\n')



