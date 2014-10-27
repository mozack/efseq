package efseq;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

public class CountTags {

	public static void run(String filename, int[] thresholds, String tag) {
		Map<Integer, Integer> buckets = new HashMap<Integer, Integer>();
		for (int threshold : thresholds) {
			buckets.put(threshold, 0);
		}
		
        File inputFile = new File(filename);
        
        SAMFileReader inputSam = new SAMFileReader(inputFile);
        inputSam.setValidationStringency(ValidationStringency.SILENT);

        for (SAMRecord read : inputSam) {
        	int value = read.getAttribute(tag) != null ? read.getIntegerAttribute(tag) : -1;
        	
        	if (value > -1) {
        		for (int threshold : thresholds) {
        			if (value <= threshold) {
        				buckets.put(threshold, buckets.get(threshold) + 1);
        				break;
        			}
        		}
        	}
        }
        
        StringBuffer out = new StringBuffer();
        for (int threshold : thresholds) {
        	out.append(threshold);
        	out.append('\t');
        }
        
        System.out.println(out.toString());
        
        out = new StringBuffer();
        for (int threshold : thresholds) {
        	out.append(buckets.get(threshold));
        	out.append('\t');
        }
        
        System.out.println(out.toString());
        
        inputSam.close();
	}
	
	public static void main(String[] args) {
		String file = args[0];
		String tag = args[1];
		String thresholdStr = args[2];
		
//		String file = "/home/lmose/dev/efseq/piotr_test1/sscs.bam";
//		String tag = "ZS";
//		String thresholdStr = "3,4,5,10,25,100,100000";
		
//		String file = "/home/lmose/dev/efseq/piotr_test1/sscs.bam";
//		String tag = "ZA";
//		String thresholdStr = "0,1,2,3,5,10,25,50,100000";
		
		
		String[] thresholdStrings = thresholdStr.split(",");
		int[] thresholds = new int[thresholdStrings.length];
		int i = 0;
		for (String threshold : thresholdStrings) {
			thresholds[i++] = Integer.parseInt(threshold);
		}		
		
		run(file, thresholds, tag);
	}
}
