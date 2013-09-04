efseq
=====
prep 
----
* Use prior to alignment to pull tag information into the read name.  
* Read pairs containing tags that lack sufficient complexity are discarded.
* Sequence and qualities longer than the max read length param are truncated

Sample usage:

java -Xmx1G -jar efseq.jar prep 1.fastq 2.fastq out1.fastq out2.fastq 100


dcs
---
* Identifies consensus sequence from a "prep"ped and aligned set of reads.
* Reads having tags with less than the specified max_variance are collapsed
* Generates intermediate sscs bam file in addition to the final dcs bam.
* Does NOT currently attempt to "pair mates" for final DCS records.
* SSCS SAM tags
	* ZS - number of reads collapsed to this SSCS record
* DCS SAM tags
	* ZD - number of SSCS collapsed to DCS (should generally be 2)
	* Z1 - original number of reads contributing to SSCS read 1
	* Z2 - original number of reads contributing to SSCS read 2

Sample usage:

java -Xmx4G -jar efseq.jar dcs sorted.bam sscs.bam dcs.bam 10