#ifndef PARSE_H
#define PARSE_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>

// struct for storing Fastq records
struct FastqRecord {
    std::string header;
    std::string sequence;
    std::string quality;
};
// struct for storing sequence records
struct SeqRecord {
    std::string sequence;
};

std::pair<std::vector<FastqRecord>, std::vector<SeqRecord>> parseFastq(const std::string& filename);

#endif
