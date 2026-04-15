#include "parse.cpp"
#include "algorithm.cpp"

int main() {
    std::string filename = "Bioinformatika - jeleni-2/fastq/J32_B_CE_IonXpress_008.fastq";

    auto result = parseFastq(filename);
    auto allRecords = result.first;
    auto seqRecords = result.second;

    for (int i = 0; i < std::min(3, (int)allRecords.size()); i++) {
        std::cout << allRecords[i].header << std::endl;
        std::cout << allRecords[i].sequence << " (duljina: " << allRecords[i].sequence.length() << ")" << std::endl;
        std::cout << std::endl;
    }

    Parameters params = {1, -1, -2};

    AlignmentResult alignment = needlemanWunsch(seqRecords[0].sequence, seqRecords[1].sequence, params);
    std::cout << "Alignment score: " << alignment.score << std::endl;
    std::cout << "Aligned Seq 1: " << alignment.alignedSeq1 << std::endl;
    std::cout << "Aligned Seq 2: " << alignment.alignedSeq2 << std::endl;
    
    return 0;
}