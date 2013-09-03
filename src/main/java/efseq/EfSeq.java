package main.java.efseq;

import java.io.IOException;

public class EfSeq {

	public static void run(String[] args) throws IOException {
		if (args.length < 1) {
			usage();
		} else if (args[0].equals("prep")) {
			String fastq1 = args[1];
			String fastq2 = args[2];
			String out1 = args[3];
			String out2 = args[4];
			int maxReadLen = Integer.parseInt(args[5]);
			
			FastqPreprocessor fp = new FastqPreprocessor();
			fp.preProcess(fastq1, fastq2, out1, out2, maxReadLen);
		} else {
			usage();
		}
	}
	
	public static void main(String[] args) throws Exception {
		run(args);
	}
	
	private static void usage() {
		System.out.println("foo");
	}
}
