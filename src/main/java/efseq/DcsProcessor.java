package efseq;

import java.io.File;
import java.io.IOException;
import java.util.List;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

/**
 * Entry point for Duplex Consensus Sequence (DCS) processing
 * 
 * @author Lisle Mose
 */
public class DcsProcessor {

	private int PHRED_SCALED_LOD_PLACEHOLDER = 70;
	
	public void dcs(String inputBam, String sscsOutput, String dcsOutput, int maxVarianceSscs, int maxVarianceDcs,
			int maxReadsAtLocus) throws IOException {

        SamLocusReader reader = new SamLocusReader(inputBam, maxReadsAtLocus);
        
        SAMFileHeader header = reader.getFileHeader();
        header.setSortOrder(SortOrder.coordinate);
        
		SAMFileWriter sscsWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				header, true, new File(sscsOutput));
		
		SAMFileWriter dcsWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				header, true, new File(dcsOutput));

		int sscsCount = 0;
		int dcsCount = 0;
		
        for (List<SAMRecord> reads : reader) {
        	
        	ConsensusSequence cs = new ConsensusSequence();
        	List<SAMRecord> sscsReads = cs.collapse(reads, maxVarianceSscs, 3, false, PHRED_SCALED_LOD_PLACEHOLDER);
        	
        	for (SAMRecord sscsRead : sscsReads) {
        		sscsWriter.addAlignment(sscsRead);
        	}
        	
        	sscsCount += sscsReads.size();
        	
        	List<SAMRecord> dcsReads = cs.collapse(sscsReads, maxVarianceDcs, 2, true, PHRED_SCALED_LOD_PLACEHOLDER);
        	
        	for (SAMRecord dcsRead : dcsReads) {
        		dcsWriter.addAlignment(dcsRead);
        	}
        	
        	dcsCount += dcsReads.size();
        }
        
        System.out.println("Closing output BAMs");
        
        sscsWriter.close();
        dcsWriter.close();
        reader.close();
        
        System.out.println("Num SSCS reads: " + sscsCount);
        System.out.println("Num DCS reads: " + dcsCount);
	}
	
	public static void main(String[] args) throws Exception {
//		String in = "/home/lmose/dev/dcs/temp1_all.sort.bam";
//		String out = "/home/lmose/dev/dcs/dcs_candidates.bam";
		
//		String in = "/home/lmose/dev/dcs/efseq/pre.chr1.bam";
//		String sscs = "/home/lmose/dev/dcs/efseq/sscs2.bam";
//		String dcs = "/home/lmose/dev/dcs/efseq/dcs2.bam";
		
		String in = "/home/lmose/dev/dcs/efseq/chr17.bam";
		String sscs = "/home/lmose/dev/dcs/efseq/sscs17.bam";
		String dcs = "/home/lmose/dev/dcs/efseq/dcs17.bam";

		
		DcsProcessor dcsProcessor = new DcsProcessor();
		dcsProcessor.dcs(in, sscs, dcs, 0, 10, 1000);
	}
}
