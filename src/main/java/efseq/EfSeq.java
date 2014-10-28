package efseq;

import java.io.IOException;

/**
 * Main entry point for efseq.
 * 
 * @author Lisle Mose
 */
public class EfSeq {

	public static void run(String[] args) throws IOException {
		if (args.length < 1) {
			usage();
		} else if (args[0].equals("prep") && args.length == 6) {
			String fastq1 = args[1];
			String fastq2 = args[2];
			String out1 = args[3];
			String out2 = args[4];
			int maxReadLen = Integer.parseInt(args[5]);
			
			FastqPreprocessor fp = new FastqPreprocessor();
			fp.preProcess(fastq1, fastq2, out1, out2, maxReadLen);
		} else if (args[0].equals("sscs") && args.length == 7) {
			String input = args[1];
			String sscs = args[2];
			int maxVarianceSscs = Integer.parseInt(args[3]);
			int maxReadsAtLocus = Integer.parseInt(args[4]);
			int phredScaledLod = Integer.parseInt(args[5]);
			int minNumReads = Integer.parseInt(args[6]);
			
			SscsProcessor sscsProcessor = new SscsProcessor();
			sscsProcessor.sscs(input, sscs, maxVarianceSscs, maxReadsAtLocus, phredScaledLod, minNumReads);			
		} else if (args[0].equals("dcs") && args.length == 7) {
			// String inputBam, String sscsOutput, String dcsOutput, int maxVariance
			String input = args[1];
			String sscs = args[2];
			String dcs = args[3];
			int maxVarianceSscs = Integer.parseInt(args[4]);
			int maxVarianceDcs = Integer.parseInt(args[5]);
			int maxReadsAtLocus = Integer.parseInt(args[6]);
			
			DcsProcessor dcsProcessor = new DcsProcessor();
			dcsProcessor.dcs(input, sscs, dcs, maxVarianceSscs, maxVarianceDcs, maxReadsAtLocus);
		} else {
			usage();
		}
	}
	
	public static void main(String[] args) throws Exception {
		run(args);
		
//		String[] argz = {
//				"sscs", "/home/lmose/dev/efseq/piotr_test1/brca2.bam", "/home/lmose/dev/efseq/piotr_test1/sscs_brca2.bam", "2", "10000", "70", "5"
//		};
//		run(argz);
	}
	
	private static void usage() {
		System.out.println("prep = Fastq Prep for alignment and downstream DCS.  sscs = Single Strand Consensus Sequence.  dcs = Duplex Consensus Sequence");
		System.out.println("java -jar efseq.jar prep <1.fastq> <2.fastq> <out1.fastq> <out2.fastq> <max_read_len>");
		System.out.println("java -jar efseq.jar sscs <input.bam> <sscs.bam> <max_variance_sscs> <max_reads_at_locus> <min_phred_scaled_lod>");
		System.out.println("java -jar efseq.jar dcs <input.bam> <sscs.bam> <dcs.bam> <max_variance_sscs> <max_variance_dcs> <max_reads_at_locus>");
	}
}
