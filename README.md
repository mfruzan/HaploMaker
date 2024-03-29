# HaploMaker

## Installation

HaploMaker is cross platform and does not require installation, it only needs java 1.8 or higher environment. Just create a directory and copy MFbio.jar into it.

Files required for pipeline:
1) Fasta reference file and its index (For example arabidopsis.fa and arabidopsis.fa.fai)
2) sorted bam file (For example Chi.bam)
3) VCF file matches bam file (For example chi.vcf)

```bash
java -Xmx20g -jar MFbio.jar --task diploidhap --vcf chi.vcf --out haplotypes.hap --afl 400 --seqtype pairedend --minmapq 10 --maxmapq 50 --ref arabidopsis.fa --refx arabidopsis.fa.fai --bam Chi.bam > out.log
```

## Parameter Description

`vcf`: VCF file</br>
`out`: output haplotype file</br>
`ref`: Reference file</br>
`refx`: Reference file index</br> 
`bam`: sorted BAM file</br>
`afl`: Average Fragment Length in base pair (In paired-end reads use average of fragments' size in the library, When using single end reads, just use read length here, when using PacBio reads use average length of reads)</br>
`seqtype`: Use `pairedend` for illumina paired-end reads, `clr` for pacbio CLR or subreads and `hifi` for Pacbio HiFi or CCS reads</br>
`minmapq`: Minimum read quality MAPQ, reads having quality less than `minmapq` will be ignored.</br>
`maxmapq`: Maximum MapQ , it depends on aligner, for Blasr use 50, for BWA, Minimap2 and pbmm2 use 60, for Bowtie2 use 42 </br>
`Xmx` is Java related parameter and controls maximum amount of memory used by Java, adjust it based on memory installed on the machine.


