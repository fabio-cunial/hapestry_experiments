import java.io.*;


/**
 * Given a dipcall VCF, the program:
 * - adds SVTYPE and SVLEN;
 * - removes short calls;
 * - removes non-PASS calls and sets FILTER to PASS.
 *
 * Remark: the program assumes that the VCF is biallelic, which is not always
 * true in dipcall. To make a VCF biallelic, use $bcftools norm$.
 */
public class FilterDipcall {
    /**
     * 
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_VCF = args[0];
        final int MIN_SV_LENGTH = Integer.parseInt(args[1]);
        final String OUTPUT_VCF = args[2];
        
        int i;
        int refLength, altLength, delta;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        String[] tokens;
        
        br = new BufferedReader(new FileReader(INPUT_VCF));
        bw = new BufferedWriter(new FileWriter(OUTPUT_VCF));
        str=br.readLine();
        while (str!=null) {
            if (str.charAt(0)=='#') {
                if (str.substring(0,6).equalsIgnoreCase("#CHROM")) {
                    bw.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variation\">\n");
                    bw.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variation\">\n");
                }
                bw.write(str); bw.newLine();
                str=br.readLine();
                continue;
            }
            tokens=str.split("\t");
            if (!tokens[6].equalsIgnoreCase(".") && !tokens[6].equalsIgnoreCase("PASS")) {
                str=br.readLine();
                continue;
            }
            tokens[6]="PASS";
            refLength=tokens[3].length(); altLength=tokens[4].length();
            delta=altLength-refLength;
            if (delta<0) delta=-delta;
            if (delta<MIN_SV_LENGTH) {
                str=br.readLine();
                continue;
            }
            if (altLength>refLength) tokens[7]="SVTYPE=INS;SVLEN="+delta;
            else if (altLength<refLength) tokens[7]="SVTYPE=DEL;SVLEN="+delta;
            bw.write(tokens[0]);
            for (i=1; i<tokens.length; i++) bw.write("\t"+tokens[i]);
            bw.newLine();
            str=br.readLine();
        }
        br.close(); bw.close();
    }
    
}