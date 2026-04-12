import java.util.Arrays;
import java.io.*;


/**
 * Given two BED files with format `CHR,FIRST,LAST,DISTANCE` computed from 
 * vcfdist's `superclusters.tsv` (each BED containing disjoint intervals), the 
 * program builds maximal connected components of overlapping or adjacent 
 * intervals across the two files, and for each component it prints the sum of
 * all distances in each file.
 */
public class CompareSuperclusters {
    /**
     * Prints to STDOUT one record per component, and to STDERR the fraction of 
     * components whose left distance is smaller that the right distance.
     */
	public static void main(String[] args) throws IOException {
		final String LEFT_BED = args[0];
        final int N_LEFT = Integer.parseInt(args[1]);
        final String RIGHT_BED = args[2];
        final int N_RIGHT = Integer.parseInt(args[3]);
        
        final int N_INTERVALS = N_LEFT+N_RIGHT;
		
        int i, j, k, p;
        int componentLast, sumLeft, sumRight, denominator;
        double numeratorLeft, numeratorRight;
        String str;
		BufferedReader br;
        String[] tokens;
        Supercluster[] intervals;
        
        // Loading all intervals from both files into the same array
        intervals = new Supercluster[N_INTERVALS];
        p=-1;
        br = new BufferedReader(new FileReader(LEFT_BED));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            intervals[++p] = new Supercluster(true,tokens[0],Integer.parseInt(tokens[1]),Integer.parseInt(tokens[2]),Integer.parseInt(tokens[3]));
            str=br.readLine();
        }
        br.close();
        br = new BufferedReader(new FileReader(RIGHT_BED));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            intervals[++p] = new Supercluster(false,tokens[0],Integer.parseInt(tokens[1]),Integer.parseInt(tokens[2]),Integer.parseInt(tokens[3]));
            str=br.readLine();
        }
        br.close();
        Arrays.sort(intervals,0,N_INTERVALS);
        
        // Computing connected components of overlapping or adjacent intervals
        numeratorLeft=0.0; numeratorRight=0.0; denominator=0;
        i=0;
        while (i<N_INTERVALS) {
            componentLast=intervals[i].last;
            for (j=i+1; j<N_INTERVALS; j++) {
                if (!intervals[j].chr.equals(intervals[i].chr) || intervals[j].first>componentLast) break;
                if (intervals[j].last>componentLast) componentLast=intervals[j].last;
            }
            sumLeft=0; sumRight=0;
            for (k=i; k<j; k++) {
                if (intervals[k].isLeft) sumLeft+=intervals[k].distance;
                else sumRight+=intervals[k].distance;
            }
            System.out.println(intervals[i].chr+"\t"+intervals[i].first+"\t"+componentLast+"\t"+sumLeft+"\t"+sumRight);
            denominator++;
            if (sumLeft<sumRight) numeratorLeft+=1.0;
            else if (sumRight<sumLeft) numeratorRight+=1.0;
            i=j;
        }
        System.err.println((numeratorLeft/denominator)+","+((denominator-numeratorLeft-numeratorRight)/denominator)+","+(numeratorRight/denominator));
	}
    
    
    /**
     * From the left BED or from the right BED
     */
    private static class Supercluster implements Comparable {
        public boolean isLeft;
        public int first, last, distance;
        public String chr;
        
        public Supercluster(boolean i, String c, int f, int l, int d) {
            this.isLeft=i;
            this.chr=c;
            this.first=f;
            this.last=l;
            this.distance=d;
        }
        
        public int compareTo(Object other) {
            Supercluster otherSupercluster = (Supercluster)other;
            int n = chr.compareTo(otherSupercluster.chr);
            if (n!=0) return n;
            if (first<otherSupercluster.first) return -1;
            else if (first>otherSupercluster.first) return 1;
            else return 0;
        }
    }

}