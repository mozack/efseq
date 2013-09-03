package main.java.efseq;

import java.io.File;
import java.io.IOException;
import java.util.List;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

/**
 * Entry point for Duplex Consensus Sequence (DCS) processing
 * 
 * @author Lisle Mose
 */
public class DcsProcessor {

	public void dcs(String inputBam, String sscsOutput, String dcsOutput, int maxVariance) throws IOException {

        SamLocusReader reader = new SamLocusReader(inputBam);
        
		SAMFileWriter sscsWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				reader.getFileHeader(), false, new File(sscsOutput));
		
		SAMFileWriter dcsWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				reader.getFileHeader(), false, new File(dcsOutput));

		int sscsCount = 0;
		int dcsCount = 0;
		
        for (List<SAMRecord> reads : reader) {
        	
        	ConsensusSequence cs = new ConsensusSequence();
        	List<SAMRecord> sscsReads = cs.collapse(reads, maxVariance, 3, false);
        	
        	for (SAMRecord sscsRead : sscsReads) {
        		sscsWriter.addAlignment(sscsRead);
        	}
        	
        	sscsCount += sscsReads.size();
        	
        	List<SAMRecord> dcsReads = cs.collapse(sscsReads, maxVariance, 2, true);
        	
        	for (SAMRecord dcsRead : dcsReads) {
        		dcsWriter.addAlignment(dcsRead);
        	}
        	
        	dcsCount += dcsReads.size();
        }
        
        sscsWriter.close();
        dcsWriter.close();
        reader.close();
        
        System.out.println("Num SSCS reads: " + sscsCount);
        System.out.println("Num DCS reads: " + dcsCount);
	}
	
	public static void main(String[] args) throws Exception {
//		String in = "/home/lmose/dev/dcs/temp1_all.sort.bam";
//		String out = "/home/lmose/dev/dcs/dcs_candidates.bam";
		
		String in = "/home/lmose/dev/dcs/efseq/pre.chr1.bam";
		String sscs = "/home/lmose/dev/dcs/efseq/sscs.bam";
		String dcs = "/home/lmose/dev/dcs/efseq/dcs.bam";
		
		DcsProcessor dcsProcessor = new DcsProcessor();
		dcsProcessor.dcs(in, sscs, dcs, 10);
	}
}
