#include "parse.h"

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
            if (currentRecord.sequence.length() == 296) {
                records.push_back(currentRecord);
                seqrecords.push_back(currentSeqRecord);
            }
        }
    }
    
    file.close();
    return {records, seqrecords};
}
