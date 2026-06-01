#include "parse.h"

// function that parses a Fastq file and returns a pair made up of a vector of FastqRecord 
// and a vector of SeqRecord
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
            // we're only interested in records where the sequence length is between 291 and 301
            if (currentRecord.sequence.length() <= 301 && currentRecord.sequence.length() >= 291) {
                records.push_back(currentRecord);
                seqrecords.push_back(currentSeqRecord);
            }
        }
    }
    
    file.close();
    return {records, seqrecords};
}
