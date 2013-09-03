efseq
=====
prep 
* Use prior to alignment to pull tag information into the read name.  
* Read pairs containing tags that lack sufficient complexity are discarded.
* Sequence and qualities longer than the max read length param are truncated

Sample usage:
java -Xmx1G -jar efseq.jar prep 1.fastq 2.fastq out1.fastq out2.fastq 100


dcs
* Identifies consensus sequence from a "prep"ped and aligned set of reads.
* Reads having tags with less than the specified max_variance are collapsed
* Generates intermediate sscs bam file in addition to the final dcs bam.

Sample usage:

java -Xmx1G -jar efseq.jar dcs sorted.bam sscs.bam dcs.bam 10