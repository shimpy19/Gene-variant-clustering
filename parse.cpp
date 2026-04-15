#include "parse.h"
#include <iostream>
#include <fstream>
#include <utility>


std::pair<std::vector<FastqRecord>, std::vector<SeqRecord>> parseFastq(const std::string& filename) {
    std::vector<FastqRecord> records;
    std::vector<SeqRecord> seqrecords;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Greška: Ne mogu otvoriti datoteku: " << filename << std::endl;
        return {records, seqrecords};
    }
    
    std::string line;
    int lineNum = 0;
    FastqRecord currentRecord;
    SeqRecord currentSeqRecord;
    
    while (std::getline(file, line)) {
        lineNum++;
        int position = lineNum % 4;
        
        if (position == 1) {
            // Header
            currentRecord.header = line;
        }
        else if (position == 2) {
            // Sekvenca (DNA)
            currentRecord.sequence = line;
            currentSeqRecord.sequence = line;
        }
        else if (position == 3) {
            // Separator (+)
        }
        else if (position == 0) {
            // Quality
            currentRecord.quality = line;
        }
        if (currentRecord.sequence.length() >= 290 && currentRecord.sequence.length() <= 305) {
                records.push_back(currentRecord);
                seqrecords.push_back(currentSeqRecord);
            }
    }
    
    file.close();
    return {records, seqrecords};
}

// Testni main, ne moras gledat
int main() {
    std::string filename = "Bioinformatika - jeleni-2/fastq/J32_B_CE_IonXpress_008.fastq";
    
    auto [allRecords, seqRecords] = parseFastq(filename);
    std::cout << allRecords.size() << std::endl;

    for (int i = 0; i < std::min(3, (int)allRecords.size()); i++) {
        std::cout << allRecords[i].header << std::endl;
        std::cout << allRecords[i].sequence << " (duljina: " << allRecords[i].sequence.length() << ")" << std::endl;
        std::cout << std::endl;
    }
    
    return 0;
}
