#ifndef ALGORITHM_H
#define ALGORITHM_H
#include <string>
#include <vector>
#include <algorithm>

struct Parameters {
    int matchScore;
    int mismatchPenalty;
    int gapPenalty;
};

struct AlignmentResult {
    int score;
};

AlignmentResult needlemanWunsch(
    const std::string& seq1, const std::string& seq2, const Parameters& params
);

#endif