# VCF2HAP
Construct gene haplotype with VCF format files and snpEff annotation.

## 1. Haplotype format defined in this pipeline
Haplotype created by *gatk_vcf_to_haplotype.pl* programe contains following files:
- *prefix.haplotype* containing 6 columns: **symbol** **haplotype rank** **haplotype** **annotation** **count** **samples**
- *prefix.vars* containing 7 columns: **chrom** **position** **reference allele** **alternative allele** **covered count** **alternative count** **annotation**

## 2. Haplotype format transforming to other data table
- *haplotype_to_popart.pl*: transforming haplotype to nexus format specified for [PopART](http://popart.otago.ac.nz/index.shtml), enabling elegent haplotype network construction and plotting.
