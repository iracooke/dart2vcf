# DArT To VCF Conversion

*WARNING:* This tool is under development. 

Converts DArT 2 row formatted csv files to vcf. 

### Usage

The only required argument is a path to your dart 2row csv file.  

```bash
dart2vcf.py Report_DHal_some_dart_file.csv > dartfile.vcf
```

If you provide a reference genome with the `-g` option, `dart2vcf` will use `bwa mem` to map the DArT tag sequences to the genome and then infer genomic position of each SNP. This information is then used to construct the CHROM and POS fields of the resulting vcf file. 

```bash
dart2vcf.py -g genome_ref.fasta Report_DHal_some_dart_file.csv > dartfile.vcf
```

When using this option your genome ref should already be indexed, which you would do with `bwa` like this

```bash
bwa index genome_ref.fasta
```

DArT provides SNP information and raw counts in separate files.  If you want individual read counts for each genotype call you can use the option `-c` to provide the matching `RawCounts` file.

```bash
dart2vcf.py -g genome_ref.fasta Report_DHal_some_dart_file_SNP.csv  dart2vcf.py -c Report_DHal_some_dart_file_RawCounts.csv > dartfile.vcf
```

If you have data represented co-analysis across multiple batches there is a high chance of sample name collision.  To make it easier to properly identify samples you can use the `-b` option to prepend the batch ID to each sample

```bash
dart2vcf.py -g genome_ref.fasta Report_DHal_some_dart_file_SNP.csv  dart2vcf.py -c Report_DHal_some_dart_file_RawCounts.csv -b > dartfile.vcf
```

### Working with the vcf

The vcf produced by this tool will contain all of the DArT locus info in the INFO field of the vcf. Genotypes are encoded as `a/b` where `a` is the allele (`1` or `0`) in the first row and `b` is the allele in the second row.  

The vcf output from this tool is sorted alphabetically by chromosome/scaffold name and by position (within each chrom). 

Chromosome/scaffold names are not added to the vcf header.  You can fix this issue by indexing the vcf

```bash
bgzip mydart.vcf
tabix mydart.vcf.gz
bcftools view -h mydart.vcf
```

### Requirements

[bwa](https://github.com/lh3/bwa) is required for mapping genomic coordinates

The [htslib](http://www.htslib.org/download/) suite of tools is recommended for working with the resulting vcf. 


### TODO

- [ ] Infer ALT and REF alleles properly. Currently these are set to maj/min respectively but could be inferred from the `AlleleSequence` field
- [ ] Parse CIGAR string to handle cases involving indels
- [ ] Use a tmp file for tag sequences so this can be cleaned up properly after running `bwa`
