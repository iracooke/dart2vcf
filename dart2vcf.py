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

parser = argparse.ArgumentParser(description='Convert DArT 2 line csv to vcf')

parser.add_argument('incsv',metavar='FILE',nargs='?',type=argparse.FileType('r'),help='DaRT csv file',default=sys.stdin)

parser.add_argument('-c', '--counts',metavar='FILE',type=argparse.FileType('r'),help='RawCounts file that matches variant calls',required=False)

parser.add_argument('-o', '--vcfout',metavar='FILE',type=argparse.FileType('w'),help='Filename for vcf output',default=sys.stdout)

parser.add_argument('-g', '--genome',type=str,help='Genome fasta. Used to find genomic coordinates',required=False)

parser.add_argument('-b', '--use-batch',action='store_true',help='Disambiguate samples by prepending batch ids to sample IDs')

parser.add_argument('-m', '--ignore-pos',action='store_true',help='Ignore SNP position. Just place the SNP in the middle of the mapped region.')

args  = parser.parse_args()


# First read the header
n_skip=0
samples_start_col=0
info_names = None
batch_names = None
batch_re = re.compile("^D.*-[0-9]+")

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
		if batch_re.search(line_values[samples_start_col]):
			batch_names=line_values[samples_start_col:]
	else:
		info_names=line_values
		break;

log.info("Skipped "+str(n_skip)+" header lines. Samples start at column "+str(samples_start_col))

# If provided, read the header of the counts file and check that it is the same as the allele file
if args.counts:
	log.info("Reading RawCounts")
	n_skip_check=0
	for cline in args.counts:
		line_values = cline.strip().split(",")

		if line_values[0]=="*":
			n_skip_check+=1
		else:
			# Check that info and sample names are the same in both files
			if info_names[0:samples_start_col]!=line_values[0:samples_start_col]:
				log.error("Column names are different in genotype and RawCounts files")
				print("Columns in 2row file",info_names)
				print("Columns in counts file",line_values)
				exit()
			if n_skip_check!=n_skip:

				log.warn("Number of header lines are different in genotype and RawCounts files")
#				exit()
			break;



# convert column headers to a dictionary

info_dict = {cn: ci for ci,cn in enumerate(info_names[0:samples_start_col])}


# What type of data is this. DArTSeq of DArTTag

data_type = "DTS" if info_dict.get('AlleleID') else "DTT"


#
# Look for column indexes of essential columns
#
# The ID is usually given by AlleleID but in DArTTag data it is given by MarkerName
snp_id_c = info_dict['AlleleID'] if info_dict.get('AlleleID') is not None else info_dict['MarkerName']

# Used for genome mapping
# The trimmed sequence is preferred here since adaptors will lead to clipping or mapping failure
#
allele_sequence_c = info_dict['TrimmedSequence'] if info_dict.get('TrimmedSequence') is not None else info_dict['AlleleSequence']

# Only present in some file types
#
snp_pos_c = info_dict.get('SnpPosition',None)

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



raw_sample_names = [sn.replace(" ","").replace("\"","") for sn in info_names[samples_start_col:]]

if args.use_batch:
	if batch_names is None:
		log.error("Use batch option specified but I was unable to find row for batch names")
		exit()
	raw_sample_names = [ bn+"_"+sn for bn,sn in zip(batch_names,raw_sample_names)]


sample_names = []
sample_tally = []
for s in raw_sample_names:
	suff=""
	ndups=sample_tally.count(s)
	if ndups>0:
		suff="_"+str(ndups)
	sample_names.append(s+suff)
	sample_tally.append(s)




n_samples=len(info_names[samples_start_col:])

vcf_records={}

allele_re = re.compile("([ATGC])>([ATGC])")


rawcounts={}
if args.counts:
	for line1,line2 in itertools.zip_longest(*[args.counts]*2):
		line1_values = line1.strip().split(",")
		line2_values = line2.strip().split(",")
		ID=line1_values[snp_id_c]
		rawcounts[ID]=[line1_values[samples_start_col:],line2_values[samples_start_col:]]

for line1,line2 in itertools.zip_longest(*[args.incsv]*2):
	line1_values = line1.strip().split(",")
	line2_values = line2.strip().split(",")

	ref_seq = line1_values[allele_sequence_c]
	snp_seq = line2_values[allele_sequence_c]

	ID=line1_values[snp_id_c]


	m = allele_re.search(ID)
	if m is None:
		log.warning("Unable to parse alleles for "+ID+" skipping this SNP")
		continue;
	allele_maj = m.group(1)
	allele_min = m.group(2)

	# This gets snp positions relative to given allele sequences in a 0-based coordinate system. 
	# We keep this for now because it simplifies book keeping checks
	# 
	snp_positions = [i for i, (ref, snp) in enumerate(zip(ref_seq,snp_seq)) if ((ref != snp) and len(set([ref,snp]) & set([allele_maj,allele_min]))==2 )  ]

	if len(snp_positions)==0:
		log.warning("Ref and Snp AlleleSequences are identical for "+ID+" skipping this SNP")
		continue;

	snp_position=-1

	if args.ignore_pos:
		snp_position=int(len(ref_seq)/2)
	else:
		if snp_pos_c is not None:
			reported_snp_pos = int(line2_values[snp_pos_c])
			if len(set(snp_positions) & set([reported_snp_pos]))!=1:
				log.error("Reported SNP position "+str(reported_snp_pos)+" does not match actual snp positions in AlleleSequence at "+str(snp_positions)+" skipping this SNP")
				continue;

			snp_position=reported_snp_pos
		elif (len(snp_positions) > 1 ):
			log.warning("Multiple potential SNP positions for "+ID+" only the first is recorded in the vcf")
			snp_position = snp_positions[0]
		else:
			snp_position = snp_positions[0]


	if snp_position==-1:
		if args.genome:
			log.error("Unable to determine SNP position")
			exit()
		else:
			snp_position=1

	# Record position and add 1 to go from 0-based to 1-based (vcf) coords
	POS=str(snp_position + 1)


	# Placeholder values for genomic coordinates. If -g is specified these will be updated later
	#
	# In case they are not replaced (-g not specified) they need cleaning up
	#
	# Contig names follow the same rules as the SAM format’s reference sequence names: they may contain any printable
	# ASCII characters in the range [!-~] apart from ‘\ , "‘’ () [] {} <>’ and may not start with ‘*’ or ‘=’. 
	# Since DArT loci names typically contain <> we need to replace those
	#
	CHROM=ID.replace(">",".gt.").replace("<",".lt.")


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



	geno_1=line1_values[samples_start_col:]
	geno_2=line2_values[samples_start_col:]

	genotypes=[]
	FORMAT="GT:AD"

	if args.counts:
		if rawcounts.get(ID)==None:
			log.error("No counts for "+ID)
			exit()
		counts_a,counts_b=rawcounts[ID]
		genotypes = [ a1.replace('-','.')+"/"+a2.replace('-','.')+":"+str(c1)+","+str(c2) for a1,a2,c1,c2 in zip(geno_1,geno_2,counts_a,counts_b) ]		
	else:
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

# DArT IDs look like this 26525701|F|0-13:T>A-13
# **Mostly/Sometimes/Usually/Maybe** This indicates that the SNP position is at 13 in 0-based system
# BUT 
# ... sometimes the real position is reported in a separate column and it's different
# ... sometimes the SNP position refers to an original sequence that is not included in the file so it is meaningless
# ... sometimes the position is just plain wrong (Yes seriously!).
# SO
# Honestly all you really have is comparing the allele and ref sequences in the file. 
# BUT
# ... often there are multiple SNPs per tag. So in that case what are the alleles? Not the single bp major and minor allele reported.
# SO
# In conclusion you can
# 1. Get the raw data and do the calls yourself
# 2. Don't get too excited about precise positions. Near enough is good enough is OK right. That is the DArT way. 
# 


	for vcfid,vcfrow in vcf_records.items():
		mappings = mapping_records.get(vcfid,[])
		if len(mappings)==0:
			n_unmapped+=1
		elif len(mappings)>1:
			n_dups+=1
		elif not args.genome and perfect_cigar_re.search(mappings[0][5]) is None:
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
				base_pos = int(mappings[0][3]) 

				# Here we need the snp pos in 0-based coords to subtract 1 to convert back
				snp_pos = int(vcfrow[1]) - 1

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
		str(n_indels)+" skipped due to complex CIGAR strings")
else:
	retained_vcfrows=vcf_records.values()

# Writing the vcf
vcfout=args.vcfout


# Construct the vcf header

vcfout.write('##fileformat=VCFv4.1'+'\n')
vcfout.write('##source=dart2vcf'+'\n')
vcfout.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate Total Depth">'+'\n')

for infokey in optional_info_indexes.keys():
	info_type,info_multiplicity,info_desc = optional_info[infokey]
	vcfout.write('##INFO=<ID='+infokey+",Number="+str(info_multiplicity)+",Type="+info_type+",Description="+info_desc+'>'+'\n')


vcfout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'+'\n')

if args.counts:
	vcfout.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele count for the ref and alt alleles in the order listed">'+'\n')
else:
	vcfout.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed. Only locus average is provided">'+'\n')	

vcf_samples_row = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
vcf_samples_row.extend(sample_names)

vcfout.write('\t'.join(vcf_samples_row)+'\n')

s1=sorted(retained_vcfrows,key=lambda row: int(row[1]))
s2=sorted(s1, key=itemgetter(0))


for vcfrow in s2:
	vcfout.write('\t'.join(vcfrow)+'\n')



