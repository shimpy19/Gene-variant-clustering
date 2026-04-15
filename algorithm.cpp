#include "algorithm.h"

AlignmentResult needlemanWunsch(
    const std::string& seq1, const std::string& seq2, const Parameters& params)
{
    int n = seq1.size();
    int m = seq2.size();

    std::vector<std::vector<int>> matrix(n+1, std::vector<int>(m+1, 0));

    // initialization of the first row and column
    for(int i = 1; i <= n; i++) {
        matrix[i][0] = i * params.gapPenalty;
    }
    for(int j = 1; j <= m; j++){
        matrix[0][j] = j * params.gapPenalty;
    }

    // filling the matrix
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= m; j++){
            int matchMismatchScore;
            if(seq1[i-1] == seq2[j-1]){
                matchMismatchScore = matrix[i-1][j-1] + params.matchScore;
            } else {
                matchMismatchScore = matrix[i-1][j-1] + params.mismatchPenalty;
            }
            int gapSeq1Score = matrix[i-1][j] + params.gapPenalty;
            int gapSeq2Score = matrix[i][j-1] + params.gapPenalty;
            matrix[i][j] = std::max({matchMismatchScore, gapSeq1Score, gapSeq2Score});
        }
    }

    // retracing the optimal aliggnment
    std::string alignedSeq1, alignedSeq2;
    int i = n, j = m;
    while(i > 0 || j > 0){
        if(i > 0 && j > 0 && 
            matrix[i][j] == matrix[i-1][j-1] + 
            (seq1[i-1] == seq2[j-1] ? params.matchScore : params.mismatchPenalty)){
            alignedSeq1 = seq1[i-1] + alignedSeq1;
            alignedSeq2 = seq2[j-1] + alignedSeq2;
            i--;
            j--;
        } else if(i > 0 && matrix[i][j] == matrix[i-1][j] + params.gapPenalty){
            alignedSeq1 = seq1[i-1] + alignedSeq1;
            alignedSeq2 = '-' + alignedSeq2;
            i--;
        } else {
            alignedSeq1 = '-' + alignedSeq1;
            alignedSeq2 = seq2[j-1] + alignedSeq2;
            j--;
        }
    }

    return {matrix[n][m], alignedSeq1, alignedSeq2};
}