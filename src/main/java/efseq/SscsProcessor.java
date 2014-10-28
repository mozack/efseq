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
 * Entry point for single strand consensus sequence processor
 * 
 * @author Lisle Mose
 */
public class SscsProcessor {

	public void sscs(String inputBam, String sscsOutput, int maxVarianceSscs,
			int maxReadsAtLocus, int phredScaledLod, int minNumReads) throws IOException {

        SamLocusReader reader = new SamLocusReader(inputBam, maxReadsAtLocus);
        
        SAMFileHeader header = reader.getFileHeader();
        header.setSortOrder(SortOrder.coordinate);
        
		SAMFileWriter sscsWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				header, true, new File(sscsOutput));

		int sscsCount = 0;
		
        for (List<SAMRecord> reads : reader) {
        	
        	ConsensusSequence cs = new ConsensusSequence();
        	List<SAMRecord> sscsReads = cs.collapse(reads, maxVarianceSscs, minNumReads, false, phredScaledLod);
        	
        	for (SAMRecord sscsRead : sscsReads) {
        		sscsWriter.addAlignment(sscsRead);
        	}
        	
        	sscsCount += sscsReads.size();
        }
        
        System.out.println("Closing output BAM");
        
        sscsWriter.close();
        reader.close();
        
        System.out.println("Num SSCS reads: " + sscsCount);
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

		
		SscsProcessor sscsProcessor = new SscsProcessor();
//		dcsProcessor.dcs(in, sscs, dcs, 0, 10, 1000);
	}
}
