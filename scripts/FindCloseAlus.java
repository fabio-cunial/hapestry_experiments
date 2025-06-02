import java.util.*;
import java.io.*;
import java.util.zip.*;


/**
 * 
 */
public class FindCloseAlus {
    
    private static final int CAPACITY = 1000;
    private static final int LENGTH_THRESHOLD = 50;
    private static int MIN_DISTANCE, SIZEMAX;
    
    /**
     * 
     *
     * @param args 
     */
    public static void main(String[] args) throws IOException {
        final String VCF_GZ = args[0];
        MIN_DISTANCE=Integer.parseInt(args[1]);  // 100
        SIZEMAX=Integer.parseInt(args[2]);  // 10000
        
        boolean found;
        int i, j, p, q;
        int currentFirst, last, pos, pos1, pos2, length1, length2, delta;
        String str, chr, currentChr;
        BufferedReader br;
        String[] window = new String[CAPACITY];
        
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(VCF_GZ))));
        str=br.readLine();
        while (str.charAt(0)=='#') str=br.readLine();
        last=-1; currentChr=""; currentFirst=-1;
        while (str!=null) {
            p=str.indexOf("\t"); q=str.indexOf("\t",p+1);
            chr=str.substring(0,p);
            pos=Integer.parseInt(str.substring(p+1,q));
            if (!chr.equals(currentChr) || pos>currentFirst+SIZEMAX) {
                processWindow(window,last);
                currentChr=chr; currentFirst=pos; last=0;
            }
            else {
                last++;
                if (last==window.length) {
                    String[] newWindow = new String[window.length<<1];
                    System.arraycopy(window,0,newWindow,0,window.length);
                    window=newWindow;
                }
            }
            window[last]=str;
            str=br.readLine();
        }
        br.close();
        processWindow(window,last);
    }
    
    
    private static void processWindow(String[] window, int last) {
        boolean found, occurs1, occurs2;
        int i, j, p, q;
        int pos1, pos2, length1, length2, delta, from, to;
        String chr;
        String[] tokens;
        
        found=false;
        for (i=0; i<last; i++) {
            p=window[i].indexOf("\t"); q=window[i].indexOf("\t",p+1);
            pos1=Integer.parseInt(window[i].substring(p+1,q));
            for (j=i+1; j<=last; j++) {
                p=window[j].indexOf("\t"); q=window[j].indexOf("\t",p+1);
                pos2=Integer.parseInt(window[j].substring(p+1,q));
                if (pos2>=pos1+MIN_DISTANCE) { 
                    tokens=window[i].split("\t"); length1=tokens[4].length(); occurs1=occurs(tokens);
                    tokens=window[j].split("\t"); length2=tokens[4].length(); occurs2=occurs(tokens);
                    delta=length1-length2;
                    if (delta<0) delta=-delta;
                    if (delta<=LENGTH_THRESHOLD && occurs1 && occurs2) { found=true; break; }
                }
            }
            if (found) break;
        }
        if (found) {
            tokens=window[0].split("\t"); from=Integer.parseInt(tokens[1]);
            tokens=window[last].split("\t"); to=Integer.parseInt(tokens[1]); chr=tokens[0];
            System.out.println(chr+"\t"+from+"\t"+to);
            // System.out.println("Window with "+(last+1)+" Alus:");
            // for (i=0; i<=last; i++) System.out.println(window[i]);
            // System.out.println();
        }
    }
    
    
    private static final boolean occurs(String[] tokens) {
        int i, p;
        
        for (i=9; i<tokens.length; i++) {
            if (tokens[i].charAt(0)=='1' || tokens[i].charAt(2)=='1') return true;
        }
        return false;
    }
    
}