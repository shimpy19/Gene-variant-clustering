#ifndef ALGORITHM_H
#define ALGORITHM_H
#include <string>
#include <vector>
#include <algorithm>

// struct for alignment parameters
struct Parameters {
    int matchScore;
    int mismatchPenalty;
    int gapPenalty;
};

// struct for alignment result
struct AlignmentResult {
    int score;
};

AlignmentResult needlemanWunsch(
    const std::string& seq1, const std::string& seq2, const Parameters& params
);

#endif