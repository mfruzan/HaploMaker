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
`afl`: DNA fragment length (In paired-end reads use average of fragments' size in the library, When using single end reads, just use read length here, when using PacBio reads use average length of reads)</br>
`seqtype`: 'pairedend' for illumina paired-end reads, 'clr' for pacbio CLR or subreads, 'hifi' for Pacbio HiFi or CCS reads</br>
`minmapq`: Minimum read quality MAPQ, reads having quality less than minmapq will be ignored.</br>
`maxmapq`: Maximum MapQ , it depends on aligner, for Blasr use 50 and for Minimap2 use 60</br>


