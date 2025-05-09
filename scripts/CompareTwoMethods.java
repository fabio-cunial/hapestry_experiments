import java.util.*;
import java.io.*;


/**
 * 
 */
public class CompareTwoMethods {
    
    /**
     *
     * @param args 
     */
    public static void main(String[] args) throws IOException {
        final String ROOT_DIR = args[0];
        final double MIN_VALUE = Double.parseDouble(args[1]);
        final double MAX_VALUE = Double.parseDouble(args[2]);
        final double BIN_LENGTH = Double.parseDouble(args[3]);
        final String TOOL_1 = args[4];
        final String TOOL_2 = args[5];
        final int MEAN_TYPE = Integer.parseInt(args[6]);  // 0=arithmetic, 1=geometric
        
        final String FILE_ID = "nodes.csv";  // File with input values
        final int TOKEN_ID = 4;  // Column with input values
        
        final int N_BINS = (int)Math.ceil((MAX_VALUE-MIN_VALUE)/BIN_LENGTH);
        int i, j, x, y;
        int denominator1, denominator2;
        long oneBetterThanTwo, twoBetterThanOne, total;
        double numerator1, numerator2;
        String str;
        File subdir;
        BufferedReader br;
        String[] directories, tokens;
        long[][] matrix;
        
        // Building the matrix
        matrix = new long[N_BINS][N_BINS];
        directories=(new File(ROOT_DIR)).list();
        for (i=0; i<directories.length; i++) {
            subdir = new File(ROOT_DIR+"/"+directories[i]);
            if (!subdir.isDirectory()) continue;
            // Tool 1
            br = new BufferedReader(new FileReader(ROOT_DIR+"/"+directories[i]+"/"+TOOL_1+"/"+FILE_ID));
            str=br.readLine(); str=br.readLine();  // Skipping header
            numerator1=MEAN_TYPE==0?0.0:1.0; denominator1=0;
            while (str!=null) {
                tokens=str.split(",");
                if (MEAN_TYPE==0) numerator1+=Double.parseDouble(tokens[TOKEN_ID]); 
                else if (MEAN_TYPE==1) numerator1*=Double.parseDouble(tokens[TOKEN_ID]); 
                denominator1++;
                str=br.readLine();
            }
            br.close();
            // Tool 2
            br = new BufferedReader(new FileReader(ROOT_DIR+"/"+directories[i]+"/"+TOOL_2+"/"+FILE_ID));
            str=br.readLine(); str=br.readLine();  // Skipping header
            numerator2=MEAN_TYPE==0?0.0:1.0; denominator2=0;
            while (str!=null) {
                tokens=str.split(",");
                if (MEAN_TYPE==0) numerator2+=Double.parseDouble(tokens[TOKEN_ID]);
                else if (MEAN_TYPE==1) numerator2*=Double.parseDouble(tokens[TOKEN_ID]);
                denominator2++;
                str=br.readLine();
            }
            br.close();
            // Updating the matrix
            x=-1; y=-1;
            if (MEAN_TYPE==0) {
                x=(int)(((numerator1/denominator1)-MIN_VALUE)/BIN_LENGTH);
                if (x>=N_BINS) x=N_BINS-1;
                if (x<0) x=0;
                y=(int)(((numerator2/denominator2)-MIN_VALUE)/BIN_LENGTH);
                if (y>=N_BINS) y=N_BINS-1;
                if (y<0) y=0;
            }
            else if (MEAN_TYPE==1) {
                x=(int)((Math.pow(numerator1,1.0/denominator1)-MIN_VALUE)/BIN_LENGTH);
                if (x>=N_BINS) x=N_BINS-1;
                if (x<0) x=0;
                y=(int)((Math.pow(numerator2,1.0/denominator2)-MIN_VALUE)/BIN_LENGTH);
                if (y>=N_BINS) y=N_BINS-1;
                if (y<0) y=0;
            }
            matrix[x][y]++;
        }
        
        // Simple stats
        total=0; oneBetterThanTwo=0; twoBetterThanOne=0;
        for (i=0; i<N_BINS; i++) {
            for (j=0; j<i; j++) { oneBetterThanTwo+=matrix[i][j]; total+=matrix[i][j]; }
            total+=matrix[i][i];
            for (j=i+1; j<N_BINS; j++) { twoBetterThanOne+=matrix[i][j]; total+=matrix[i][j]; }
        }
        System.err.println(TOOL_1+" is better than "+TOOL_2+" in "+oneBetterThanTwo+" windows ("+((100.0*oneBetterThanTwo)/total)+"%)");
        System.err.println(TOOL_2+" is better than "+TOOL_1+" in "+twoBetterThanOne+" windows ("+((100.0*twoBetterThanOne)/total)+"%)");
        System.err.println("Total windows: "+total);
        
        // Printing the matrix
        for (i=0; i<N_BINS; i++) {
            for (j=0; j<N_BINS; j++) System.out.print(matrix[i][j]+",");
            System.out.println();
        }
    }
    
}