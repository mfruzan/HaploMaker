# HaploMiner
Files required for pipeline:
1) Fasta reference file and its index (For example arabidopsis.fa and arabidopsis.fa.fai)
2) sorted bam file (For example Chi.bam)


java -Xmx3000m -jar MFbio.jar --task=diploidhap --vcf=chi.vcf --out haplotypes.hap --p1 400 --p2 2 --p3 10 --p4 2 --ref arabidopsis.fa --refx arabidopsis.fa.fai --bam Chi.bam > out.log </br>
Parameter Description:
vcf: VCF file
out: output haplotype file
ref: Reference file
refx: Reference file index 
bam: sorted BAM file
p1: Maximum insert size (When using single end reads, just read length is used here)
p2: 1 for single end, 2 for paired end reads
p3: Minimum read quality MAPQ
p4: minimum haplotype size (in number of SNPs) to be reported


