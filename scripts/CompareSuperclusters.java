import java.util.Arrays;
import java.io.*;


/**
 * Given two BED files with format `CHR,FIRST,LAST,DISTANCE` extracted from 
 * vcfdist's `superclusters.tsv` (each BED contains disjoint intervals), the 
 * program builds maximal connected components of overlapping or adjacent 
 * intervals across the two files, and for each component it prints the sum of
 * all distances in each file.
 */
public class CompareSuperclusters {
    
    /**
     *
     */
	public static void main(String[] args) throws IOException {
		final String LEFT_BED = args[0];
        final int N_LEFT = Integer.parseInt(args[1]);
        final String RIGHT_BED = args[2];
        final int N_RIGHT = Integer.parseInt(args[3]);
        
        final int N_INTERVALS = N_LEFT+N_RIGHT;
		
        int i, j, k, p;
        int clusterLast, sumLeft, sumRight;
        String str;
		BufferedReader br;
        String[] tokens;
        Interval[] intervals;
        
        // Loading all intervals from both files into a single list
        intervals = new Interval[N_INTERVALS];
        p=-1;
        br = new BufferedReader(new FileReader(LEFT_BED));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            intervals[++p] = new Interval(true,tokens[0],Integer.parseInt(tokens[1]),Integer.parseInt(tokens[2]),Integer.parseInt(tokens[3]));
            str=br.readLine();
        }
        br.close();
        br = new BufferedReader(new FileReader(RIGHT_BED));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            intervals[++p] = new Interval(false,tokens[0],Integer.parseInt(tokens[1]),Integer.parseInt(tokens[2]),Integer.parseInt(tokens[3]));
            str=br.readLine();
        }
        br.close();
        Arrays.sort(intervals,0,N_INTERVALS);
        
        // Computing connected components of overlapping or adjacent intervals
        i=0;
        while (i<N_INTERVALS) {
            clusterLast=intervals[i].last;
            for (j=i+1; j<N_INTERVALS; j++) {
                if (!intervals[j].chr.equals(intervals[i].chr) || intervals[j].first>clusterLast) break;
                if (intervals[j].last>clusterLast) clusterLast=intervals[j].last;
            }
            sumLeft=0; sumRight=0;
            for (k=i; k<j; k++) {
                if (intervals[k].isLeft) sumLeft+=intervals[k].distance;
                else sumRight+=intervals[k].distance;
            }
            System.out.println(intervals[i].chr+"\t"+intervals[i].first+"\t"+clusterLast+"\t"+sumLeft+"\t"+sumRight);
            i=j;
        }
	}
    
    
    private static class Interval implements Comparable {
        public boolean isLeft;
        public int first, last, distance;
        public String chr;
        
        public Interval(boolean i, String c, int f, int l, int d) {
            this.isLeft=i;
            this.chr=c;
            this.first=f;
            this.last=l;
            this.distance=d;
        }
        
        public int compareTo(Object other) {
            Interval otherInterval = (Interval)other;
            int n = chr.compareTo(otherInterval.chr);
            if (n!=0) return n;
            if (first<otherInterval.first) return -1;
            else if (first>otherInterval.first) return 1;
            else return 0;
        }
    }

}