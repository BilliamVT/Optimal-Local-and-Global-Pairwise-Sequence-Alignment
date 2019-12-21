import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;

public class Main {
    // Main method
    public static void main(String args[]) {
        // Initialize sequence1 and sequence 2
        String sequence1 = "";
        String sequence2 = "";

        // Initialize input file
        File inputFile = new File("2.in");

        // Read in from file
        try {
            // Scanner to read from file
            Scanner fileIn = new Scanner(inputFile);

            // Set sequence1 and sequence2 to the sequences in the file
            sequence1 = fileIn.nextLine();
            sequence2 = fileIn.nextLine();

            // Close fileIn
            fileIn.close();

            // If file not found
        } catch (FileNotFoundException e) {
            // Print file not found
            System.out.println("Input file not found");
        }

        // Initialize match, misMatch, gapPenalty
        int match = 2;
        int misMatch = -1;
        int gapPenalty = -2;

        // Initialize optimalScore
        int optimalScore = -99;

        // Initialize optimalAlignment
        String optimalAlignment[] = new String[]{};

        // Initialize multipleOptimal (0 for no 1 for yes)
        boolean multipleOptimal = false;

        // Initialize alignmentMatrix
        int[][] alignmentMatrix = {};

        // User chooses whether to run global or local alignment
        Scanner userInput = new Scanner(System.in);
        System.out.println("Would you like to run Needleman-Wunsch (0) or Smith-Waterman (1) pairwise alignment: ");
        int alignment = userInput.nextInt();
        userInput.close();

        if (alignment == 0) {
            System.out.println("Running global alignment...");
            // Returns matrix after global alignment
            alignmentMatrix = NeedlemanWunsch(sequence1, sequence2, match, misMatch, gapPenalty);

            // Get optimal score from bottom right corner
            optimalScore = alignmentMatrix[sequence1.length()][sequence2.length()];

            // Get if multiple optimal
            multipleOptimal = FindIfMultipleOptimalGlobal(alignmentMatrix, sequence1, sequence2, match, misMatch, gapPenalty);

            // Find optimal alignment for global
            optimalAlignment = FindGlobalOptimal(alignmentMatrix, sequence1, sequence2, match, misMatch, gapPenalty);

        } else if (alignment == 1) {
            System.out.println("Running local alignment...");

            // Returns matrix after local alignment
            alignmentMatrix = SmithWaterman(sequence1, sequence2, match, misMatch, gapPenalty);

            // Get optimal score from bottom right corner
            optimalScore = alignmentMatrix[sequence1.length()][sequence2.length()];

            // Get if multiple optimal
            multipleOptimal = FindIfMultipleOptimalLocal(alignmentMatrix, sequence1, sequence2, match, misMatch, gapPenalty);

            // Find optimal alignment for local
            optimalAlignment = FindLocalOptimal(alignmentMatrix, sequence1, sequence2, match, misMatch, gapPenalty);

        }

        /***************************************************************************************************************
         ***************************************************************************************************************
         * Begin write to files
         ***************************************************************************************************************
         ***************************************************************************************************************/
        // Write optimal score to file
        try {
            // Create writer
            BufferedWriter writer1 = new BufferedWriter(new FileWriter("2.o1"));

            // Write optimal score to file
            writer1.write(Integer.toString(optimalScore));

            // Close writer
            writer1.close();

            // Catch IOException
        } catch (IOException e) {
            System.out.println("IOException error writing to file 2.o1");
        }

        // Write dynamic programming matrix to file
        try {
            // Create writer
            BufferedWriter writer2 = new BufferedWriter(new FileWriter("2.o2"));

            // Print alignment matrix to file
            for (int i = 0; i < sequence1.length() + 1; i++) {
                // So there isn't a blank row at the top
                if (i != 0) {
                    writer2.write("\n");
                }
                // Write each node
                for (int j = 0; j < sequence2.length() + 1; j++) {
                    writer2.write(alignmentMatrix[i][j] + " ");
                }
            }

            // Close writer
            writer2.close();

            // Catch IOException
        } catch (IOException e) {
            System.out.println("IOException error writing to file 2.o2");
        }

        // Write one optimal alignment to file
        try {
            // Create writer
            BufferedWriter writer3 = new BufferedWriter(new FileWriter("2.o3"));

            // Write optimal alignment to file
            writer3.write(optimalAlignment[0]);
            writer3.write("\n");
            writer3.write(optimalAlignment[1]);

            // Close writer
            writer3.close();

            // Catch IOException
        } catch (IOException e) {
            System.out.println("IOException error writing to file 2.o3");
        }

        // Write YES or NO depending on if there are multiple optimal alignments
        try {
            // Create writer
            BufferedWriter writer4 = new BufferedWriter(new FileWriter("2.o4"));

            // Write YES or NO
            if (multipleOptimal) {
                writer4.write("YES");
            } else {
                writer4.write("NO");
            }


            // Close writer
            writer4.close();

            // Catch IOException
        } catch (IOException e) {
            System.out.println("IOException error writing to file 2.o4");
        }
    }

    /*******************************************************************************************************************
     *******************************************************************************************************************
     * End write to files
     *******************************************************************************************************************
     *******************************************************************************************************************/

    /*******************************************************************************************************************
     *******************************************************************************************************************
     * Needleman-Wunsch Global Alignment
     *******************************************************************************************************************
     *******************************************************************************************************************/
    // Runs Needleman-Wunsch global alignment algorithm
    static int[][] NeedlemanWunsch(String sequence1, String sequence2, int match, int misMatch, int gapPenalty) {
        // Each sequence starts with 0
        sequence1 = "0" + sequence1;
        sequence2 = "0" + sequence2;

        // Make matrix length of sequence1 by length of sequence2
        int[][] matrix = new int[sequence1.length()][sequence2.length()];

        // Add gap penalty along edges
        for (int i = 0; i < sequence1.length(); i++) {
            matrix[i][0] = i * gapPenalty;
        }
        for (int i = 0; i < sequence2.length(); i++) {
            matrix[0][i] = i * gapPenalty;
        }

        // Run alignment
        for (int i = 1; i < sequence1.length(); i++) {
            for (int j = 1; j < sequence2.length(); j++) {
                // Find score for gap on way
                int gapOne = matrix[i - 1][j] + gapPenalty;

                // Find score for gap other way
                int gapTwo = matrix[i][j - 1] + gapPenalty;

                // Initialize diagonal
                int diagonal;

                // Find score for diagonal
                if (sequence1.charAt(i) == sequence2.charAt(j)) {
                    diagonal = matrix[i - 1][j - 1] + match;
                } else {
                    diagonal = matrix[i - 1][j - 1] + misMatch;
                }

                // Find max of three scores
                int max = Math.max(Math.max(gapOne, gapTwo), diagonal);

                // Set score for point i,j
                matrix[i][j] = max;
            }
        }

        // Return filled in matrix
        return matrix;
    }

    /*******************************************************************************************************************
     *******************************************************************************************************************
     * Smith-Waterman Local Alignment
     *******************************************************************************************************************
     *******************************************************************************************************************/
    // Runs Smith-Waterman local alignment algorithm
    static int[][] SmithWaterman(String sequence1, String sequence2, int match, int misMatch, int gapPenalty) {
        // Each sequence starts with 0
        sequence1 = "0" + sequence1;
        sequence2 = "0" + sequence2;

        // Make matrix length of sequence1 by length of sequence2
        int[][] matrix = new int[sequence1.length()][sequence2.length()];

        // Add 0s along edges
        for (int i = 0; i < sequence1.length(); i++) {
            matrix[i][0] = 0;
        }
        for (int i = 0; i < sequence2.length(); i++) {
            matrix[0][i] = 0;
        }

        // Run alignment
        for (int i = 1; i < sequence1.length(); i++) {
            for (int j = 1; j < sequence2.length(); j++) {
                // Find score for gap on way
                int gapOne = matrix[i - 1][j] + gapPenalty;

                // Find score for gap other way
                int gapTwo = matrix[i][j - 1] + gapPenalty;

                // Initialize diagonal
                int diagonal;

                // Find score for diagonal
                if (sequence1.charAt(i) == sequence2.charAt(j)) {
                    diagonal = matrix[i - 1][j - 1] + match;
                } else {
                    diagonal = matrix[i - 1][j - 1] + misMatch;
                }

                // Find max of three scores and zero because of local alignment
                int max = Math.max(Math.max(Math.max(gapOne, gapTwo), diagonal), 0);

                // Set score for point i,j
                matrix[i][j] = max;
            }
        }

        // Return filled in matrix
        return matrix;
    }

    /*******************************************************************************************************************
     *******************************************************************************************************************
     * Global Optimal Alignment
     *******************************************************************************************************************
     *******************************************************************************************************************/
    // Finds optimal sequence for global alignment
    static String[] FindGlobalOptimal(int[][] matrix, String sequence1, String sequence2, int match, int misMatch, int gapPenalty) {
        // Sequences start with 0
        sequence1 = "0" + sequence1;
        sequence2 = "0" + sequence2;

        // Stringbuilder of the alignments
        StringBuilder align1 = new StringBuilder();
        StringBuilder align2 = new StringBuilder();

        // Starting point
        int i = sequence1.length() - 1;
        int j = sequence2.length() - 1;

        // While not at the top left corner
        while (i != 0 && j != 0) {
            // Score based off of match or misMatch
            int score;

            // If they match set score to match
            if (sequence1.charAt(i) == sequence2.charAt(j)) {
                score = match;
            }
            // Else set to misMatch
            else {
                score = misMatch;
            }

            // If score equals score coming from gap
            if (matrix[i][j] == (matrix[i][j - 1]) + gapPenalty) {
                // Insert gap
                align1.insert(0, "-");
                align2.insert(0, sequence2.charAt(j));

                // Only decrement j as it only goes in that direction
                j--;

                // If score equals coming from diagonal
            } else if (matrix[i][j] == (matrix[i - 1][j - 1] + score)) {
                // Insert both amino acids
                align1.insert(0, sequence1.charAt(i));
                align2.insert(0, sequence2.charAt(j));

                // Decrement both i and j to go diagonal
                i--;
                j--;

                // If score equals coming from other gap
            } else if (matrix[i][j] == (matrix[i - 1][j] + gapPenalty)) {
                // Insert gap
                align1.insert(0, sequence1.charAt(i));
                align2.insert(0, "-");

                // Only decrement i as it only goes in that direction
                i--;
            }
        }

        // Returning both sequences
        return new String[]{align1.toString(), align2.toString()};
    }

    /*******************************************************************************************************************
     *******************************************************************************************************************
     * Local Optimal Alignment
     *******************************************************************************************************************
     *******************************************************************************************************************/
    // Finds optimal sequence for local alignment
    static String[] FindLocalOptimal(int[][] matrix, String sequence1, String sequence2, int match, int misMatch, int gapPenalty) {
        // Must find start node
        // Initialize max score
        int maxScore = 0;

        // Must record starting location of max score
        int startingI = 0;
        int startingJ = 0;

        // Loop through to find starting location and maxScore
        for (int i = 0; i < sequence1.length() + 1; i++) {
            for (int j = 0; j < sequence2.length() + 1; j++) {
                if (matrix[i][j] >= maxScore) {
                    maxScore = matrix[i][j];
                    startingI = i;
                    startingJ = j;
                }
            }
        }

        // Sequences start with 0
        sequence1 = "0" + sequence1;
        sequence2 = "0" + sequence2;

        // Stringbuilder of the alignments
        StringBuilder align1 = new StringBuilder();
        StringBuilder align2 = new StringBuilder();

        // Starting point
        int i = startingI;
        int j = startingJ;

        // While not at the top left corner or at a zero
        while (i != 0 && j != 0 && matrix[i][j] != 0) {
            // Score based off of match or misMatch
            int score;

            // If they match set score to match
            if (sequence1.charAt(i) == sequence2.charAt(j)) {
                score = match;
            }
            // Else set to misMatch
            else {
                score = misMatch;
            }

            // If score equals score coming from gap
            if (matrix[i][j] == (matrix[i][j - 1]) + gapPenalty) {
                // Insert gap
                align1.insert(0, "-");
                align2.insert(0, sequence2.charAt(j));

                // Only decrement j as it only goes in that direction
                j--;

                // If score equals coming from diagonal
            } else if (matrix[i][j] == (matrix[i - 1][j - 1] + score)) {
                // Insert both amino acids
                align1.insert(0, sequence1.charAt(i));
                align2.insert(0, sequence2.charAt(j));

                // Decrement both i and j to go diagonal
                i--;
                j--;

                // If score equals coming from other gap
            } else if (matrix[i][j] == (matrix[i - 1][j] + gapPenalty)) {
                // Insert gap
                align1.insert(0, sequence1.charAt(i));
                align2.insert(0, "-");

                // Only decrement i as it only goes in that direction
                i--;
            }
        }

        // Returning both sequences
        return new String[]{align1.toString(), align2.toString()};
    }

    /*******************************************************************************************************************
     *******************************************************************************************************************
     * Finds if there are multiple optimal alignments for global alignment
     *******************************************************************************************************************
     *******************************************************************************************************************/
    static boolean FindIfMultipleOptimalGlobal(int[][] matrix, String sequence1, String sequence2, int match, int misMatch, int gapPenalty){
        // Sequences start with 0
        sequence1 = "0" + sequence1;
        sequence2 = "0" + sequence2;

        // Stringbuilder of the alignments
        StringBuilder align1 = new StringBuilder();
        StringBuilder align2 = new StringBuilder();

        // Starting point
        int i = sequence1.length() - 1;
        int j = sequence2.length() - 1;

        // While not at the top left corner
        while (i != 0 && j != 0) {
            // Keep track of paths
            boolean[] paths = new boolean[]{false, false, false};

            // Score based off of match or misMatch
            int score;

            // If they match set score to match
            if (sequence1.charAt(i) == sequence2.charAt(j)) {
                score = match;
            }
            // Else set to misMatch
            else {
                score = misMatch;
            }

            // If score equals score coming from gap
            if (matrix[i][j] == (matrix[i][j - 1]) + gapPenalty) {
                // Path used
                paths[0] = true;

                // Insert gap
                align1.insert(0, "-");
                align2.insert(0, sequence2.charAt(j));

                // Only decrement j as it only goes in that direction
                j--;

                // If score equals coming from diagonal
            }
            if (matrix[i][j] == (matrix[i - 1][j - 1] + score)) {
                // Path used
                paths[1] = true;

                // Insert both amino acids
                align1.insert(0, sequence1.charAt(i));
                align2.insert(0, sequence2.charAt(j));

                // Decrement both i and j to go diagonal
                i--;
                j--;

                // If score equals coming from other gap
            }
            if (matrix[i][j] == (matrix[i - 1][j] + gapPenalty)) {
                // Path used
                paths[2] = true;

                // Insert gap
                align1.insert(0, sequence1.charAt(i));
                align2.insert(0, "-");

                // Only decrement i as it only goes in that direction
                i--;
            }
            if (paths[0] && paths[1] || paths[0] && paths[2] || paths[1] && paths[2]){
                return true;
            }
        }
        return false;
    }

    /*******************************************************************************************************************
     *******************************************************************************************************************
     * Finds if there are multiple optimal alignments for local alignment
     *******************************************************************************************************************
     *******************************************************************************************************************/
    static boolean FindIfMultipleOptimalLocal(int[][] matrix, String sequence1, String sequence2, int match, int misMatch, int gapPenalty){
        // Must find start node
        // Initialize max score
        int maxScore = 0;

        // Must record starting location of max score
        int startingI = 0;
        int startingJ = 0;

        // Loop through to find starting location and maxScore
        for (int i = 0; i < sequence1.length() + 1; i++) {
            for (int j = 0; j < sequence2.length() + 1; j++) {
                if (matrix[i][j] >= maxScore) {
                    maxScore = matrix[i][j];
                    startingI = i;
                    startingJ = j;
                }
            }
        }

        // Sequences start with 0
        sequence1 = "0" + sequence1;
        sequence2 = "0" + sequence2;

        // Stringbuilder of the alignments
        StringBuilder align1 = new StringBuilder();
        StringBuilder align2 = new StringBuilder();

        // Starting point
        int i = startingI;
        int j = startingJ;

        // While not at the top left corner or at a zero
        while (i != 0 && j != 0 && matrix[i][j] != 0) {
            // Keep track of paths
            boolean[] paths = new boolean[]{false, false, false};

            // Score based off of match or misMatch
            int score;

            // If they match set score to match
            if (sequence1.charAt(i) == sequence2.charAt(j)) {
                score = match;
            }
            // Else set to misMatch
            else {
                score = misMatch;
            }

            // If score equals score coming from gap
            if (matrix[i][j] == (matrix[i][j - 1]) + gapPenalty) {
                // Insert gap
                align1.insert(0, "-");
                align2.insert(0, sequence2.charAt(j));

                // Path used
                paths[0] = true;

                // Only decrement j as it only goes in that direction
                j--;

                // If score equals coming from diagonal
            }
            if (matrix[i][j] == (matrix[i - 1][j - 1] + score)) {
                // Insert both amino acids
                align1.insert(0, sequence1.charAt(i));
                align2.insert(0, sequence2.charAt(j));

                // Path used
                paths[1] = true;

                // Decrement both i and j to go diagonal
                i--;
                j--;

                // If score equals coming from other gap
            }
            if (matrix[i][j] == (matrix[i - 1][j] + gapPenalty)) {
                // Insert gap
                align1.insert(0, sequence1.charAt(i));
                align2.insert(0, "-");

                // Path used
                paths[2] = true;

                // Only decrement i as it only goes in that direction
                i--;
            }
            if (paths[0] && paths[1] || paths[0] && paths[2] || paths[1] && paths[2]){
                return true;
            }
        }
        return false;
    }

}