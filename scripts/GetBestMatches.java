import java.util.HashMap;
import java.io.*;


/**
 * Given a BED file (computed by `bedtools intersect`) that contains every 
 * possible overlap between a left set of superclusters and a right set of 
 * superclusters, the program prints every left supercluster with just a right
 * supercluster with maximum overlap.
 */
public class GetBestMatches {

	public static void main(String[] args) throws IOException {
		final String INPUT_TSV = args[0];
		
        int overlap, bestOverlap;
        String str, leftInterval, currentLeftInterval, bestMatch;
		BufferedReader br;
        String[] tokens;
        
        br = new BufferedReader(new FileReader(INPUT_TSV));
        currentLeftInterval=""; bestOverlap=0; bestMatch="";
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            leftInterval=tokens[0]+","+tokens[1]+","+tokens[2];
            overlap=Integer.parseInt(tokens[8]);
            if (leftInterval.equals(currentLeftInterval)) {
                if (overlap>bestOverlap) {
                    bestOverlap=overlap;
                    bestMatch=str;
                }
            }
            else {
                if (currentLeftInterval.length()>0) System.out.println(bestMatch);
                currentLeftInterval=leftInterval;
                bestOverlap=overlap;
                bestMatch=str;
            }
            str=br.readLine();
        }
        System.out.println(bestMatch);
        br.close();
	}

}