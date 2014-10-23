package efseq;

import java.io.IOException;

/**
 * Prepares paired end fastq files for initial alignment and subsequent DCS processing.
 * 
 * @author Lisle Mose
 */
public class FastqPreprocessor {
	
	public static final int TAG_LENGTH =  8;
	public static final int CONNECTOR_LENGTH = 13;
	public static final int SEP_LENGTH = 1;
	
	// 8 + 13 + 1
	public static final int VALID_SEQUENCE_START_IDX = 22;

	/**
	 * Prep input FASTQ files for initial alignment and downstream consensus sequence processing. 
	 */
	public void preProcess(String fastqInput1, String fastqInput2, 
			String fastqOutput1, String fastqOutput2, int maxReadLength) throws IOException {
		FastqInputFile fastq1 = new FastqInputFile();
		FastqInputFile fastq2 = new FastqInputFile();
		
		FastqOutputFile out1 = new FastqOutputFile();
		FastqOutputFile out2 = new FastqOutputFile();
		
		fastq1.init(fastqInput1);
		fastq2.init(fastqInput2);
		
		out1.init(fastqOutput1);
		out2.init(fastqOutput2);
		
		FastqRecord rec1 = fastq1.getNextRecord();
		FastqRecord rec2 = fastq2.getNextRecord();
		
		int cnt = 1;
		
		while ((rec1 != null) && (rec2 != null)) {

//			String tag = rec1.getSequence().substring(0, TAG_LENGTH) +
//					rec2.getSequence().substring(0, TAG_LENGTH);
			
			String tag1 = rec1.getSequence().substring(0, TAG_LENGTH);
			String tag2 = rec2.getSequence().substring(0, TAG_LENGTH);
			
			String connector1 = rec1.getSequence().substring(TAG_LENGTH, TAG_LENGTH+CONNECTOR_LENGTH);
			String connector2 = rec2.getSequence().substring(TAG_LENGTH, TAG_LENGTH+CONNECTOR_LENGTH);
			
			String sep1 = rec1.getSequence().substring(TAG_LENGTH+CONNECTOR_LENGTH, TAG_LENGTH+CONNECTOR_LENGTH+SEP_LENGTH);
			String sep2 = rec2.getSequence().substring(TAG_LENGTH+CONNECTOR_LENGTH, TAG_LENGTH+CONNECTOR_LENGTH+SEP_LENGTH);

			rec1.setId("@" + tag1 + "_" + connector1 + "_" + sep1 + "_" + rec1.getId());
			rec2.setId("@" + tag2 + "_" + connector2 + "_" + sep2 + "_" + rec2.getId());
			trimSequenceAndQuality(rec1, maxReadLength);
			trimSequenceAndQuality(rec2, maxReadLength);
			
			out1.write(rec1);
			out2.write(rec2);

			
			/*
			if (isTagSufficientlyComplex(tag)) {
				rec1.setId("@" + tag + "_" + rec1.getId());
				rec2.setId("@" + tag + "_" + rec2.getId());
				trimSequenceAn dQuality(rec1, maxReadLength);
				trimSequenceAndQuality(rec2, maxReadLength);
				
				out1.write(rec1);
				out2.write(rec2);
			}
			*/
			
			rec1 = fastq1.getNextRecord();
			rec2 = fastq2.getNextRecord();
			
			if ((cnt++ % 1000000) == 0) {
				System.err.println("Processed " + (cnt-1) + " read pairs.");
			}
		}
		
		out1.close();
		out2.close();
		fastq1.close();
		fastq2.close();
		
		System.err.println("Done.");
	}
	
	// Skip initial identification bases and trim off end bases which may be of low quality depending
	// upon sequencing technology.
	private void trimSequenceAndQuality(FastqRecord rec, int maxReadLength) {
		int endIdx = maxReadLength + VALID_SEQUENCE_START_IDX;
		rec.setSequence(rec.getSequence().substring(VALID_SEQUENCE_START_IDX, endIdx < rec.getSequence().length() ? endIdx : rec.getSequence().length()));
		rec.setQuality(rec.getQuality().substring(VALID_SEQUENCE_START_IDX, endIdx < rec.getQuality().length() ? endIdx : rec.getQuality().length()));		
	}
	
	// Require tag to not contain 10+ repeats of a base or an ambiguous base
	private boolean isTagSufficientlyComplex(String tag) {
		return !tag.matches(".*N.*|.*A{10}.*|.*T{10}.*|.*C{10}.*|.*G{10}.*");
	}
	
	public static void main(String[] args) {
		FastqPreprocessor fp = new FastqPreprocessor();
		System.out.println(fp.isTagSufficientlyComplex("ATCGATCGATCGATCG"));
		System.out.println(fp.isTagSufficientlyComplex("ATCNGATCGATCGATCG"));
		System.out.println(fp.isTagSufficientlyComplex("ATCGATCGGGGGGGGGGATCGATCG"));
		System.out.println(fp.isTagSufficientlyComplex("NATCGATCGATCGATCG"));
		System.out.println(fp.isTagSufficientlyComplex("ATCGATCGATCGATCGN"));
		System.out.println(fp.isTagSufficientlyComplex("ATCGATCGGGGGGGGGGGG"));
		System.out.println(fp.isTagSufficientlyComplex("AAAAAAAAAAATCGATCGATCGATCG"));
	}
}
