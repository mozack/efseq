package main.java.efseq;

import java.io.IOException;

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
		} else if (args[0].equals("dcs") && args.length == 5) {
			// String inputBam, String sscsOutput, String dcsOutput, int maxVariance
			String input = args[1];
			String sscs = args[2];
			String dcs = args[3];
			int maxVariance = Integer.parseInt(args[4]);
			
			DcsProcessor dcsProcessor = new DcsProcessor();
			dcsProcessor.dcs(input, sscs, dcs, maxVariance);
		} else {
			usage();
		}
	}
	
	public static void main(String[] args) throws Exception {
		run(args);
	}
	
	private static void usage() {
		System.out.println("prep = Fastq Prep for alignment and downstream DCS.  dcs = Duplex Consensus Sequence");
		System.out.println("java -jar efseq.jar prep <1.fastq> <2.fastq> <out1.fastq> <out2.fastq> <max_read_len>");
		System.out.println("java -jar efseq.jar dcs <input.bam> <sscs.bam> <dcs.bam> <max_variance>");
	}
}
