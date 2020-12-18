# HaploMiner

Installation:
Haplominer is cross platform and does not require installation, it only needs java 1.8 or higher environment. Just create a directory and copy MFbio.jar into it.
Files required for pipeline:
1) Fasta reference file and its index (For example arabidopsis.fa and arabidopsis.fa.fai)
2) sorted bam file (For example Chi.bam)
3) VCF file matches bam file (For example chi.vcf)


java -Xmx3000m -jar MFbio.jar --task diploidhap --vcf=chi.vcf --out haplotypes.hap --p1 400 --p2 2 --p3 10 --p4 2 --ref arabidopsis.fa --refx arabidopsis.fa.fai --bam Chi.bam > out.log </br>
Parameter Description:</br>
vcf: VCF file</br>
out: output haplotype file</br>
ref: Reference file</br>
refx: Reference file index</br> 
bam: sorted BAM file</br>
p1: Maximum insert size (In paired-end reads just use average of fragments' size in the library, When using single end reads, just use read length here)</br>
p2: 1 for single end, 2 for paired end reads</br>
p3: Minimum read quality MAPQ</br>
p4: minimum haplotype size (in number of SNPs) to be reported (This parameter only affects output, minimum value is 2)</br>


