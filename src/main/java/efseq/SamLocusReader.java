package efseq;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

/**
 * Iterates over a Sam or Bam file returning a List of SAMRecords for each starting locus.
 * Input must be sorted by coordinate.
 * 
 * @author Lisle Mose
 */
public class SamLocusReader implements Iterable<List<SAMRecord>> {
	
	private SAMFileReader inputSam;
	private Iterator<SAMRecord> iter;
	private SAMRecord cachedRead;
	private int numRecords = 0;
	private int maxReadsAtLocus;

    public SamLocusReader(String filename, int maxReadsAtLocus) {
        File inputFile = new File(filename);
        
        this.maxReadsAtLocus = maxReadsAtLocus;
        
        inputSam = new SAMFileReader(inputFile);
        inputSam.setValidationStringency(ValidationStringency.SILENT);
  
        iter = inputSam.iterator();
    }
    
    public SAMFileHeader getFileHeader() {
    	return inputSam.getFileHeader();
    }
    
    private boolean isReadUnmapped(SAMRecord read) {
    	return read.getReadUnmappedFlag() && read.getAlignmentStart() == 0;
    }
    
    private boolean hasMoreReads() {
    	if ((cachedRead == null) && (iter.hasNext())) {
    		cachedRead = advanceSamIter();
    	}
    	
    	return cachedRead != null && !isReadUnmapped(cachedRead);
    }
    
    private SAMRecord advanceSamIter() {
    	SAMRecord read = iter.next();
		
    	numRecords++;
		
		if ((numRecords % 1000000) == 0) {
			System.out.println("Processed " + numRecords + " records.  Chromosome: " + read.getReferenceName());
		}
		
    	return read;
    }
    
    private SAMRecord getNextRead() {
    	SAMRecord read = null;
    	
    	if (cachedRead != null) {
    		read = cachedRead;
    		cachedRead = null;
    	} else {
    		read = advanceSamIter();;
    	}
    	
    	return read;
    }
    
    private String getLocus(SAMRecord read) {
    	return isReadUnmapped(read) ? "unmapped" : read.getReferenceName() + ":" + read.getAlignmentStart();
    }
    
    private List<SAMRecord> getNextReadList() {
    	List<SAMRecord> reads = new ArrayList<SAMRecord>();
    	
    	SAMRecord read = null;
    	String locus = null;
    	
    	if (hasMoreReads()) {
    		read = getNextRead();
    		locus = getLocus(read);
    		    		
    		reads.add(read);
    		
    		// getLocus() must be called prior to hasMoreReads() in this check.
    		// hasMoreReads() would otherwise advance the cachedRead
    		while (getLocus(read).equals(locus) && hasMoreReads()) {
    			read = getNextRead();
    			
    			if (getLocus(read).equals(locus)) {
    				
    				if (reads.size() < maxReadsAtLocus) {
    					reads.add(read);
    				}
    			}
    		}
    		
    		// If the last read doesn't have the same locus, cache it
    		if (!getLocus(read).equals(locus)) {
    			cachedRead = read;
    		}
    	}
    	
    	if (reads.size() >= maxReadsAtLocus) {
    		SAMRecord firstRead = reads.get(0);
    		System.out.println("Too many reads at: " + firstRead.getReferenceName() + ":" + firstRead.getAlignmentStart());
    		reads.clear();
    	}
    	
    	return reads;
    }
    
    @Override
	public Iterator<List<SAMRecord>> iterator() {
		return new SamLocusIterator(this);
	}
    
    public void close() {
    	inputSam.close();
    }

	private static class SamLocusIterator implements Iterator<List<SAMRecord>> {

        private SamLocusReader reader;
        
        SamLocusIterator(SamLocusReader reader) {
            this.reader = reader;
        }
        
        @Override
        public boolean hasNext() {
        	return reader.hasMoreReads();        	
        }
        
        @Override
        public List<SAMRecord> next() {
        	return reader.getNextReadList();
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Remove not supported for SamMultiMappingIterator.");
        }
    }

}
