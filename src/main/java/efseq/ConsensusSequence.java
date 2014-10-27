package efseq;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

/**
 * Collapses reads into containing matching tags into consensus sequence.
 * Handles both single strand consensus sequence (SSCS)
 * and duplex consensus sequence (DCS).  Does not currently consider pairings
 * for DCS.  Tag matching allows for mismatches and gaps.
 * 
 * @author Lisle Mose
 */
public class ConsensusSequence {
	
	private static final double MIN_CONSENSUS_BASE_FREQUENCY = .9;
	private static final char AMBIGUOUS_BASE = 'N';
	private static final char MAX_QUAL = 'H';
	private static final char MIN_QUAL = '!';

	public List<SAMRecord> collapse(List<SAMRecord> readsAtLocus, int maxVariance,
			int minFreq, boolean isDcs, int phredScaledLod) throws IOException {

		List<List<SAMRecord>> collapsedReads = new ArrayList<List<SAMRecord>>();
        
        FuzzyStringMatch fuzzy = new FuzzyStringMatch();
        
        if (readsAtLocus.size() > 10) {
//        	System.out.println("Here...");
        	int a=1;
        }
        
    	for (SAMRecord read : readsAtLocus) {
			
//    		if (read.getReadUnmappedFlag() || read.getNotPrimaryAlignmentFlag() || (read.getFlags() & 0x800) == 0x800 || read.getCigar().getCigarElement(0).getOperator() == CigarOperator.S) {
    		if (read.getReadUnmappedFlag() || read.getNotPrimaryAlignmentFlag() || (read.getFlags() & 0x800) == 0x800) {
    			// Skip to the next read
    			// TODO: Adjust for soft clipped start upstream from here.
    			continue;
    		}
    		
    		boolean isCollapsed = false;
    		
			for (List<SAMRecord> collapsedRead : collapsedReads) {
							
				SAMRecord firstCollapsedRead = collapsedRead.get(0);
				
				if (isSameAlignment(read, firstCollapsedRead)) {
					
					String firstTag;
					String secondTag;
					    		
					firstTag = collapsedRead.get(0).getReadName().substring(0, FastqPreprocessor.TAG_LENGTH);
					
					if (isDcs) {
			    		secondTag = read.getReadName().substring(12,24) + 
			    				read.getReadName().substring(0,12);
					} else {
						secondTag = read.getReadName().substring(0, FastqPreprocessor.TAG_LENGTH);
					}
		    		
		    		int matchScore = fuzzy.scoreMatch(firstTag, secondTag);
		    		
		    		int diff = FastqPreprocessor.TAG_LENGTH - matchScore;
		    		
		    		if (diff <= maxVariance) {
		    			isCollapsed = true;
		    			collapsedRead.add(read);
		    			break;
		    		}
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
    			validCollapsedReads.add(collapse(reads, isDcs, phredScaledLod));
    		}
    	}
    	
    	return validCollapsedReads;
	}
	
	private boolean isSameAlignment(SAMRecord read1, SAMRecord read2) {
		boolean isSame = false;
		
		if (read1.getReferenceName().equals(read2.getReferenceName()) &&
			read1.getAlignmentStart() == read2.getAlignmentStart() &&
			read1.getCigar().equals(read2.getCigar()) && 
			read1.getReadNegativeStrandFlag() == read2.getReadNegativeStrandFlag() &&
			read1.getReadPairedFlag() == read2.getReadPairedFlag()) {
				
			// Now check mate info
			if (read1.getReadPairedFlag()) {
				if (read1.getMateUnmappedFlag() == read2.getMateUnmappedFlag()) {
					if (read1.getMateUnmappedFlag()) {
						// Mates for both reads are unmapped
						isSame = true;
					} else {
						if (read1.getMateReferenceName().equals(read2.getMateReferenceName()) &&
							read1.getAlignmentStart() == read2.getAlignmentStart() &&
							read1.getCigar().equals(read2.getCigar())) {
							
							isSame = true;
						}
					}
				}
			}
		}
		
		return isSame;
	}
	
	private SAMRecord collapse(List<SAMRecord> reads, boolean isDcs, int phredScaledLod) {
		
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
				
				//count += 1;
				
				// Count base qualities
				count += read.getBaseQualities()[i];
				
				counts.put(base, count);
			}
		}
		
		StringBuffer sequence = new StringBuffer();
		StringBuffer qualities = new StringBuffer();
		
		int ambiguousBaseCount = 0;
		
		for (int i=0; i<reads.get(0).getReadString().length(); i++) {
			char consensusBase = getConsensusBase(baseCounts.get(i), reads.size(), phredScaledLod);

			sequence.append(consensusBase);
			if (consensusBase == AMBIGUOUS_BASE) {
				qualities.append(MIN_QUAL);
				ambiguousBaseCount += 1;
			} else {
				qualities.append(MAX_QUAL);
			}
		}
		
		SAMRecord read = reads.get(0);
		read.setReadName(read.getReadName().substring(0, 24));
		read.setReadString(sequence.toString());
		read.setBaseQualityString(qualities.toString());
		
		if (isDcs) {
			int idx = 1;
			// Record all SCCS counts for this DCS read
			for (SAMRecord contributingRead : reads) {
				if (idx < 10) {
					int count = (Integer) contributingRead.getAttribute("ZS");
					String tag = "Z" + idx;
					read.setAttribute(tag, count);
					idx++;
				} else {
					break;
				}
			}
			
			read.setAttribute("ZD", reads.size());
			read.setAttribute("ZS", null);
			
		} else {
			read.clearAttributes();

			// Record SSCS count
			read.setAttribute("ZS", reads.size());
			read.setAttribute("ZA", ambiguousBaseCount);
		}
		
//		read.setReadPairedFlag(false);
		
		return read;
	}

	private char getConsensusBase(Map<Character, Integer> counts, int numReads, int phredScaledLod) {
		
		int totalQual = 0;
		char winningBase = AMBIGUOUS_BASE;
		int maxCount = -1;
		
		for (char base : counts.keySet()) {
			totalQual += counts.get(base);
			if (counts.get(base) > maxCount) {
				maxCount = counts.get(base);
				winningBase = base;
			}
		}
				
		int winningQual = counts.get(winningBase);
		int otherQual = totalQual - winningQual;

		if ((winningQual - otherQual) < phredScaledLod) {
			winningBase = AMBIGUOUS_BASE;
		}
		
		return winningBase;
	}

/*	
	// Check that # base observations is >= 90% of total observations
	private char getConsensusBase(Map<Character, Integer> counts, int numReads) {
		for (Character base : counts.keySet()) {
			if (counts.get(base) != null && counts.get(base) >= (numReads * MIN_CONSENSUS_BASE_FREQUENCY)) {
				return base;
			}
		}
		
		return AMBIGUOUS_BASE;
	}
*/
	
	public static void main(String[] args) throws Exception {
		String in = "/home/lmose/dev/dcs/temp1_all.sort.bam";
		String out = "/home/lmose/dev/dcs/dcs_candidates.bam";
		
		ConsensusSequence dcs = new ConsensusSequence();
//		dcs.identifyDcs(in, out, 10);
	}
}
