# HaploMiner
Files required for pipeline:
1) Fasta reference file and its index (For example arabidopsis.fa and arabidopsis.fa.fai)
2) sorted bam file (For example Chi.bam)

Steps of pipeline:

(A) create pileup file using samtools (set minimum mapq 10):
samtools mpileup -B -A -d 9999 -q 10 -f arabidpsis.fa -t DP -b bam_files.txt > chi2.pileup </br>
where bam_files.txt will contain fully qualified name of bam file.

(B) Create variant file from pileup file:
java -jar MFbio.jar showform=no task=variant srcdir=chi2.pileup destdir=genotype.var p1=1 p2=2 p3=0.1 </br>
where p1 is number of samples in pileup file, p2 is minimum number of reads for a non-reference allele (3 recommended) and p3 is minimum ratio for minor allele to call a position heterozygous. For example if we set p3=0.1 then ifd allele 1 has 95% and allele 2 has 5% then position is called homozygous, but if allele 1 is 85% and allele 2 is 15% then position is called heterozygous.

(C) Generate haplotype blocks:
java -jar MFbio.jar showform=no task=diploidhap2 srcdir=genotype.var destdir=haplotypes.hap p1=400 p2=2 p3=10 p4=10 1=arabidopsis.fa 2=arabidopsis.fa.fai file1=Chi.bam > out.log </br>
where srcdir is variant file generated at step B. destdir is name of output haplotype file. p1 is estimataded fragment length (insert size for paired-end reads, read length for single end reads). p2=1 means single end reads and p2=2 means paired end reads. p3 is minimum MAPQ (10 is recommended). p4 is minimum haplotype block size (in number of SNP in the block) to be reported. Parameter 1 refers to reference fasta file. Parameter 2 refers to reference fasta index. file1 referes to bam file.
