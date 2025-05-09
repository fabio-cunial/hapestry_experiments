import java.lang.Math;
import java.util.*;
import java.io.*;
import java.util.zip.*;


/**
 * Finds kanpig windows that contain >=2 records with similar kanpig score.
 */
public class FindKanpigWindows {
    
    /**
     * java FindKanpigWindows dipcall_50_merged.vcf.gz 50 500 0.9 0.97 4 2
     *
     * @param args names correspond to kanpig's command-line arguments.
     */
    public static void main(String[] args) throws IOException {
        final String VCF_GZ = args[0];
        final int SIZEMIN = Integer.parseInt(args[1]);  // 10
        final int NEIGHDIST = Integer.parseInt(args[2]);  // 500
        final double SIZESIM = Double.parseDouble(args[3]);  // default 0.9
        final double SEQSIM = Double.parseDouble(args[4]);  // 0.97
        final int KMER = Integer.parseInt(args[5]);  // default 4
        final int MINKFREQ = Integer.parseInt(args[6]);  // default 2
        
        final int CAPACITY = 1000;
        int last, pos, currentFirst, currentLast, nPairs;
        String str, currentChr;
        BufferedReader br;
        int[] lengths = new int[CAPACITY];
        HashMap<String,Integer>[] kmers = new HashMap[CAPACITY];
        String[] tokens;
        
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(VCF_GZ))));
        str=br.readLine();
        while (str.charAt(0)=='#') str=br.readLine();
        last=-1; currentChr=""; currentFirst=-1; currentLast=-1;
        while (str!=null) {
            tokens=str.split("\t");
            tokens[3].toLowerCase(); tokens[4].toLowerCase();
            pos=Integer.parseInt(tokens[1]);
            if (!tokens[0].equals(currentChr) || pos>currentLast+NEIGHDIST) {
                if (last>=0) {
                    nPairs=nCompatiblePairs(kmers,lengths,last,SIZEMIN,MINKFREQ,SIZESIM,SEQSIM);
                    if (nPairs!=0) System.out.println(currentChr+"\t"+currentFirst+"\t"+currentLast+"\t"+nPairs);
                }
                // Next window
                currentChr=tokens[0]; currentFirst=pos; currentLast=-1;
                last=-1;
            }
            last++;
            if (last==kmers.length) {
                HashMap<String,Integer>[] newKmers = new HashMap[kmers.length*2];
                System.arraycopy(kmers,0,newKmers,0,kmers.length);
                kmers=newKmers;
                int[] newLengths = new int[lengths.length*2];
                System.arraycopy(lengths,0,newLengths,0,lengths.length);
                lengths=newLengths;
            }
            if (kmers[last]==null) kmers[last] = new HashMap<String,Integer>();
            loadKmers(tokens[3].substring(1),tokens[4].substring(1),KMER,kmers[last]);
            if (tokens[3].length()==1 && tokens[4].length()>1) {  // INS
                lengths[last]=tokens[4].length()-1;
            }
            else if (tokens[4].length()==1 && tokens[3].length()>1) {  // DEL
                lengths[last]=tokens[3].length()-1;
            }
            else {  // Replacement
                lengths[last]=(tokens[4].length()-1+tokens[3].length()-1)/2;
            }
            if (pos+tokens[3].length()>currentLast) currentLast=pos+tokens[3].length();
            str=br.readLine();
        }
        br.close();
        if (last>=0) {
            nPairs=nCompatiblePairs(kmers,lengths,last,SIZEMIN,MINKFREQ,SIZESIM,SEQSIM);
            if (nPairs!=0) System.out.println(currentChr+"\t"+currentFirst+"\t"+currentLast+"\t"+nPairs);
        }
    }
    
    
    private static final void loadKmers(String ref, String alt, int k, HashMap<String,Integer> map) {
        int i;
        String key;
        
        map.clear();
        for (i=0; i<=alt.length()-k; i++) {
            key=alt.substring(i,i+k);
            if (map.containsKey(key)) map.put(key,Integer.valueOf(map.get(key).intValue()+1));
            else map.put(key,Integer.valueOf(1));
        }
        for (i=0; i<=ref.length()-k; i++) {
            key=ref.substring(i,i+k);
            if (map.containsKey(key)) map.put(key,Integer.valueOf(map.get(key).intValue()-1));
            else map.put(key,Integer.valueOf(-1));
        }
    }
    
    
    private static final int nCompatiblePairs(final HashMap<String,Integer>[] kmers, final int[] lengths, final int last, final int sizemin, final int minkfreq, final double sizesim, final double seqsim) {
        int i, j;
        int nPairs;
        
        nPairs=0;
        for (i=0; i<last; i++) {
            if (lengths[i]<sizemin) continue;
            for (j=i+1; j<=last; j++) {
                if (lengths[j]<sizemin) continue;
                if (Math.min(lengths[i],lengths[j])/Math.max(lengths[i],lengths[j])>=sizesim && kmerSimilarity(kmers[i],kmers[j],minkfreq)>=seqsim) nPairs++;
            }
        }
        return nPairs;
    }
    
    
    private static final double kmerSimilarity(HashMap<String,Integer> map1, HashMap<String,Integer> map2, int minkfreq) {
        int i;
        double out, count1, count2;
        String key;
        Iterator<String> iterator;
        
        out=1.0;
        iterator=map1.keySet().iterator();
        while (iterator.hasNext()) {
            key=iterator.next();
            count1=map1.get(key).intValue();
            if (map2.containsKey(key)) {
                count2=map2.get(key).intValue();
                if (Math.abs(count1)+Math.abs(count2)<minkfreq) continue;
                out-=Math.abs(count1-count2)/(Math.abs(count1)+Math.abs(count2));
            }
            else {
                if (Math.abs(count1)<minkfreq) continue;
                out-=1.0;
            }
        }
        iterator=map2.keySet().iterator();
        while (iterator.hasNext()) {
            key=iterator.next();
            count2=map2.get(key).intValue();
            if (Math.abs(count2)<minkfreq || map1.containsKey(key)) continue;
            out-=1.0;
        }
        return out;
    }
    
}