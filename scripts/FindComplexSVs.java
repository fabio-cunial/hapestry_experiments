import java.util.*;
import java.io.*;
import java.util.zip.*;


/**
 * A complex window is a set of overlapping or adjacent calls such that at least
 * two calls in the set occur on the same haplotype.
 */
public class FindComplexSVs {
    
    private static final int CAPACITY = 1000;
    private static int MAX_DISTANCE;
    
    /**
     * @param args 
     * 0: a phased truth VCF from dipcall;
     * 1: max distance between two calls for them to be considered adjacent.
     */
    public static void main(String[] args) throws IOException {
        final String VCF_GZ = args[0];
        MAX_DISTANCE=Integer.parseInt(args[1]);  // 50
        
        int p1, p2, p3, p4;
        int last, currentFirst, currentLast, firstPos, lastPos;
        String str, chr, currentChr;
        BufferedReader br;
        String[] window = new String[CAPACITY];
        
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(VCF_GZ))));
        str=br.readLine();
        while (str.charAt(0)=='#') str=br.readLine();
        last=-1; currentChr=""; currentFirst=-1; currentLast=-1;
        while (str!=null) {
            p1=str.indexOf("\t");  // CHROM
            p2=str.indexOf("\t",p1+1);  // POS
            p3=str.indexOf("\t",p2+1);  // ID
            p4=str.indexOf("\t",p3+1);  // REF
            chr=str.substring(0,p1);
            firstPos=getFirstPos(str,p1,p2,p3);
            if (!chr.equals(currentChr) || firstPos>currentLast+MAX_DISTANCE) {
                processWindow(window,last);
                currentChr=chr;
                currentLast=getLastPos(str,firstPos,p1,p2,p3,p4);
                last=0;
            }
            else {
                last++;
                if (last==window.length) {
                    String[] newWindow = new String[window.length<<1];
                    System.arraycopy(window,0,newWindow,0,window.length);
                    window=newWindow;
                }
                lastPos=getLastPos(str,firstPos,p1,p2,p3,p4);
                if (lastPos>currentLast) currentLast=lastPos;
            }
            window[last]=str;
            str=br.readLine();
        }
        br.close();
        processWindow(window,last);
    }
    
    
    private static void processWindow(String[] window, int last) {
        int i, j;
        int p1, p2, p3, p4, p5, p6, p7, p8, p9;
        int nPairs, first1, last1, first2, firstOfWindow, lastOfWindow;
        String chr, gt1, gt2;
        
        nPairs=0; chr=""; firstOfWindow=-1; lastOfWindow=-1;
        for (i=0; i<=last; i++) {
            p1=window[i].indexOf("\t");  // CHROM
            p2=window[i].indexOf("\t",p1+1);  // POS
            p3=window[i].indexOf("\t",p2+1);  // ID
            p4=window[i].indexOf("\t",p3+1);  // REF
            p5=window[i].indexOf("\t",p4+1);  // ALT
            p6=window[i].indexOf("\t",p5+1);  // QUAL
            p7=window[i].indexOf("\t",p6+1);  // FILTER
            p8=window[i].indexOf("\t",p7+1);  // INFO
            p9=window[i].indexOf("\t",p8+1);  // FORMAT
            gt1=window[i].substring(p9+1);
            if (chr.length()==0) chr=window[i].substring(0,p1);
            first1=getFirstPos(window[i],p1,p2,p3);
            if (firstOfWindow==-1) firstOfWindow=first1;
            last1=getLastPos(window[i],first1,p1,p2,p3,p4);
            if (last1>lastOfWindow) lastOfWindow=last1;
            for (j=i+1; j<=last; j++) {
                p1=window[j].indexOf("\t");  // CHROM
                p2=window[j].indexOf("\t",p1+1);  // POS
                p3=window[j].indexOf("\t",p2+1);  // ID
                first2=getFirstPos(window[j],p1,p2,p3);
                if (first2>last1+MAX_DISTANCE) break;
                p4=window[j].indexOf("\t",p3+1);  // REF
                p5=window[j].indexOf("\t",p4+1);  // ALT
                p6=window[j].indexOf("\t",p5+1);  // QUAL
                p7=window[j].indexOf("\t",p6+1);  // FILTER
                p8=window[j].indexOf("\t",p7+1);  // INFO
                p9=window[j].indexOf("\t",p8+1);  // FORMAT
                gt2=window[j].substring(p9+1);
                if (onSameHaplotype(gt1,gt2)) nPairs++;
            }
        }
        if (nPairs>0) System.out.println(chr+"\t"+firstOfWindow+"\t"+(lastOfWindow+1)+"\t"+nPairs);
    }
    
    
    /**
     * @return the first position that is affected by the SV (zero-based).
     * INS: the SV is assumed to start at `POS`.
     * DEL: the SV is assumed to start at `POS+1`.
     */
    private static final int getFirstPos(String vcfRecord, int p1, int p2, int p3) {
        final int POS = Integer.parseInt(vcfRecord.substring(p1+1,p2));
        final int POS_ZERO_BASED = POS-1;
        final int REF_LENGTH = p3-(p2+1);
        return REF_LENGTH==1?POS_ZERO_BASED:POS_ZERO_BASED+1;
    }
    
    
    /**
     * @param firstPos the first position affected by the SV;
     * @return the last position that is affected by the SV (zero-based).
     * INS: this is conventionally set to `POS+1`.
     */
    private static final int getLastPos(String vcfRecord, int firstPos, int p1, int p2, int p3, int p4) {
        final int REF_LENGTH = p3-(p2+1);
        final int SV_LENGTH = REF_LENGTH-1;
        return REF_LENGTH==1?firstPos+1:firstPos+SV_LENGTH-1;
    }
    
    
    private static final boolean onSameHaplotype(String gt1, String gt2) {
        return (gt1.charAt(0)=='1' && gt2.charAt(0)=='1') || (gt1.charAt(2)=='1' && gt2.charAt(2)=='1');
    }
    
}