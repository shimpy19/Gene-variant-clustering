#include "algorithm.h"

AlignmentResult needlemanWunsch(
    const std::string& seq1, const std::string& seq2, const Parameters& params)
{
    int n = seq1.size();
    int m = seq2.size();

    std::vector<int> prev(m+1, 0), curr(m+1, 0);

    for(int j = 1; j <= m; j++){
        prev[j] = j * params.gapPenalty;
    }

    for(int i = 1; i <= n; i++){
    curr[0] = i * params.gapPenalty;
    for(int j = 1; j <= m; j++){
        int matchMismatch = prev[j-1] + 
            (seq1[i-1] == seq2[j-1] ? params.matchScore : params.mismatchPenalty);
        int gapSeq1 = prev[j] + params.gapPenalty;
        int gapSeq2 = curr[j-1] + params.gapPenalty;
        curr[j] = std::max({matchMismatch, gapSeq1, gapSeq2});
    }

    std::swap(prev, curr);
}

return {prev[m]};
}