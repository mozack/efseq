package efseq;

/**
 * Compares 2 DNA sequences providing a score indicating similarity.
 * Allows mismatches and gaps.  This class is not threadsafe.
 * 
 * @author Lisle Mose
 */
public class FuzzyStringMatch {
	
	private String s1;
	private String s2;
	int[][] matrix;

	/**
	 * Returns a score indicating the similarity between the 2 input sequences. 
	 */
	public int scoreMatch(String s1, String s2) {
		this.s1 = s1;
		this.s2 = s2;
		matrix = new int[s1.length()][];
		for (int i=0; i<s1.length(); i++) {
			matrix[i] = new int[s2.length()];
		}
		
		for (int col=0; col<s2.length(); col++) {
			for (int row=0; row<s1.length(); row++) {
				matrix[row][col] = scoreEntry(row, col);				
			}
			
//			printMatrix();
		}
		
		int max = matrix[s1.length()-1][s2.length()-1];
		
		return max;
	}
	
	private int scoreEntry(int row, int col) {
		int val1 = getValue(row-1, col);
		int val2 = getValue(row, col-1);
		char ch1 = s1.charAt(row);
		char ch2 = s2.charAt(col);
		int val3 = getValue(row-1, col-1) + (ch1 == ch2 && ch1 != 'N' ? 1 : 0);
		
		int max = Math.max(val1, val2);
		max = Math.max(max, val3);
		
		return max;
	}
	
	private int getPrevMax(int row, int col) {
		int val1 = getValue(row-1, col);
		int val2 = getValue(row, col-1);
		int val3 = getValue(row-1, col-1);
		
		int max = Math.max(val1, val2);
		max = Math.max(max, val3);
		
		return max;
	}
	
	private int getValue(int row, int col) {
		int val = 0;
		if ((row >= 0) && (col >= 0)) {
			val = matrix[row][col];
		}
		
		return val;
	}
	
	public void printMatrix() {
		System.out.println(" " + s2);
		for (int row=0; row<s1.length();row++) {
			StringBuffer str = new StringBuffer();
			str.append(s1.charAt(row));
			for (int col=0; col<s2.length(); col++) {
				str.append(matrix[row][col]);
			}
			System.out.println(str.toString());
		}
	}
	
	public static void main(String[] args) {
		String s1 = "ATCGA";
		String s2 = "ATTCG";
		
		FuzzyStringMatch f = new FuzzyStringMatch();
		
		int score = f.scoreMatch(s1, s2);
		
		System.out.println("score: " + score);
		
		f.printMatrix();
	}
}
