# DArT To VCF Conversion

Converts DArT 2 row formatted csv files to vcf. 


### Mapping to genomic coordinates

If an appropriate reference genome is provided this script will use `bwa mem` to map the DArT tag sequences to the genome and then infer genomic position of each SNP. This information is then use to construct the CHROM and POS fields of the resulting vcf file. 
