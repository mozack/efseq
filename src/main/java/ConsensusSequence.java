package efseq;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.SAMRecord;

public class ConsensusSequence {
	
	private static final double MIN_CONSENSUS_BASE_FREQUENCY = .9;
	private static final char AMBIGUOUS_BASE = 'N';
	private static final char MAX_QUAL = 'H';
	private static final char MIN_QUAL = '!';

	public List<SAMRecord> collapse(List<SAMRecord> readsAtLocus, int maxVariance,
			int minFreq, boolean isDcs) throws IOException {

		List<List<SAMRecord>> collapsedReads = new ArrayList<List<SAMRecord>>();
        
        FuzzyStringMatch fuzzy = new FuzzyStringMatch();
        
        if (readsAtLocus.size() > 10) {
//        	System.out.println("Here...");
        	int a=1;
        }
        
    	for (SAMRecord read : readsAtLocus) {
			
    		boolean isCollapsed = false;
    		
			for (List<SAMRecord> collapsedRead : collapsedReads) {
							
				String firstTag;
				String secondTag;
				    		
				firstTag = collapsedRead.get(0).getReadName().substring(0, 24);
				
				if (isDcs) {
		    		secondTag = read.getReadName().substring(12,24) + 
		    				read.getReadName().substring(0,12);
				} else {
					secondTag = read.getReadName().substring(0, 24);
				}
	    		
	    		int matchScore = fuzzy.scoreMatch(firstTag, secondTag);
	    		
	    		int diff = 24 - matchScore;
	    		
	    		if (diff <= maxVariance) {
	    			isCollapsed = true;
	    			collapsedRead.add(read);
	    			break;
	    		}
			}
			
			if (!isCollapsed) {

				List<SAMRecord> collapsedRead = new ArrayList<SAMRecord>();
				collapsedRead.add(read);
				collapsedReads.add(collapsedRead);
			}
    	}
    	
    	List<SAMRecord> validCollapsedReads = new ArrayList<SAMRecord>();
    	for (List<SAMRecord> reads : collapsedReads) {
    		if (reads.size() >= minFreq) {
    			validCollapsedReads.add(collapse(reads, isDcs));
    		}
    	}
    	
    	return validCollapsedReads;
	}
	
	private SAMRecord collapse(List<SAMRecord> reads, boolean isDcs) {
		
		// Map of maps
		// Top map -> key = position
		// Inner map -> key = base, value = count
		Map<Integer, Map<Character, Integer>> baseCounts = new HashMap<Integer, Map<Character, Integer>>();

		for (SAMRecord read : reads) {
			
			for (int i=0; i<read.getReadString().length(); i++) {
				Map<Character, Integer> counts = baseCounts.get(i);
				if (counts == null) {
					counts = new HashMap<Character, Integer>();
					baseCounts.put(i, counts);
				}
				char base = read.getReadString().charAt(i);
				Integer count = counts.get(base);
				if (count == null) {
					count = new Integer(0);
				}
				
				count += 1;
				
				counts.put(base, count);
			}
		}
		
		StringBuffer sequence = new StringBuffer();
		StringBuffer qualities = new StringBuffer();
		
		for (int i=0; i<reads.get(0).getReadString().length(); i++) {
			char consensusBase = getConsensusBase(baseCounts.get(i), reads.size());

			sequence.append(consensusBase);
			if (consensusBase == AMBIGUOUS_BASE) {
				qualities.append(MIN_QUAL);
			} else {
				qualities.append(MAX_QUAL);
			}
		}
		
		SAMRecord read = reads.get(0);
		read.setReadString(sequence.toString());
		read.setBaseQualityString(qualities.toString());
		
		if (isDcs) {
			int idx = 1;
			// Record all SCCS counts for this DCS read
			for (SAMRecord contributingRead : reads) {
				if (idx < 10) {
					int count = (Integer) contributingRead.getAttribute("ZC");
					String tag = "Z" + idx;
					read.setAttribute(tag, count);
					idx++;
				} else {
					break;
				}
			}
			
			read.setAttribute("ZD", reads.size());
			
		} else {
			// Record SSCS count
			read.setAttribute("ZC", reads.size());
		}
		
		return read;
	}
	
	private char getConsensusBase(Map<Character, Integer> counts, int numReads) {
		for (Character base : counts.keySet()) {
			if (counts.get(base) != null && counts.get(base) >= (numReads * MIN_CONSENSUS_BASE_FREQUENCY)) {
				return base;
			}
		}
		
		return AMBIGUOUS_BASE;
	}
	
	public static void main(String[] args) throws Exception {
		String in = "/home/lmose/dev/dcs/temp1_all.sort.bam";
		String out = "/home/lmose/dev/dcs/dcs_candidates.bam";
		
		ConsensusSequence dcs = new ConsensusSequence();
//		dcs.identifyDcs(in, out, 10);
	}
}
