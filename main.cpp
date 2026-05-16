#include "parse.h"
#include "algorithm.h"
#include "clustering.h"
#include <iostream>
#include <fstream>

int main() {
    std::string filename = "Bioinformatika - jeleni-2/fastq/J29_B_CE_IonXpress_005.fastq";

    auto result = parseFastq(filename);
    auto allRecords = result.first;
    auto seqRecords = result.second;

    std::ofstream seqFile("sequences_J29.txt");
    if (seqFile.is_open()) {
        for (size_t i = 0; i < seqRecords.size(); i++) {
            seqFile << i << "\t" << seqRecords[i].sequence << "\n";
            }
            seqFile.close();
        } 
    else {
        std::cerr << "Error: Cannot open sequence file." << std::endl;
    }

    std::vector<std::vector<int>> distance_matrix(1700, std::vector<int>(1700));

    #pragma omp parallel for collapse(2)
    for(size_t i = 0; i < seqRecords.size(); i++) {
        for(size_t j = i+1; j < seqRecords.size(); j++) {
            auto alignmentResult = needlemanWunsch(seqRecords[i].sequence, seqRecords[j].sequence, {1, -1, -2});
            distance_matrix[i][j] = std::max(0, 296 - alignmentResult.score);
            distance_matrix[j][i] = std::max(0, 296 - alignmentResult.score);
        }
    }

    std::string baseName = filename.substr(filename.find_last_of("/\\") + 1);
    std::string sampleName = baseName.substr(0, baseName.find("_"));
    std::string outputFileName = "distance_matrix_" + sampleName + ".txt";
    std::string seqFileName = "sequences_" + sampleName + ".txt";

    std::ofstream outputFile(outputFileName);
    if(outputFile.is_open()){
        for(size_t i = 0; i < seqRecords.size(); i++) {
            for(size_t j = 0; j < seqRecords.size(); j++) {
                outputFile << distance_matrix[i][j] << " ";
            }
            outputFile << std::endl;
        }
        outputFile.close();
        
        runClustering(outputFileName, seqFileName);

    } else {
        std::cerr << "Error: Cannot open output file." << std::endl;
    }
    
    return 0;
}