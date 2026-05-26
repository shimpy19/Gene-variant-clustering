#include "parse.h"
#include "algorithm.h"
#include "clustering.h"
#include <iostream>
#include <fstream>
#include <filesystem>

int main() {
    namespace fs = std::filesystem;
    const std::string dirPath = "Bioinformatika - jeleni-2/fastq";

    const fs::path matricesDir = "matrices";
    const fs::path sequencesDir = "sequences";

    fs::create_directories(matricesDir);
    fs::create_directories(sequencesDir);

    for (const auto& entry : fs::directory_iterator(dirPath)) {
        if (!entry.is_regular_file()) continue;
        if (entry.path().extension() != ".fastq") continue;

        std::string filename = entry.path().string();
        std::cout << "Processing file: " << entry.path().filename().string() << std::endl;

        auto result = parseFastq(filename);
        auto allRecords = result.first;
        auto seqRecords = result.second;

        if (seqRecords.empty()) {
            std::cerr << "Warning: no sequences parsed from " << filename << std::endl;
            continue;
        }

        std::string baseName = entry.path().filename().string();
        std::string sampleName = baseName.substr(0, baseName.find("_"));
        fs::path seqFilePath = sequencesDir / ("sequences_" + sampleName + ".txt");
        fs::path matrixFilePath = matricesDir / ("distance_matrix_" + sampleName + ".txt");

        std::string seqFileName = seqFilePath.string();
        std::string outputFileName = matrixFilePath.string();

        std::ofstream seqFile(seqFileName);
        if (seqFile.is_open()) {
            for (size_t i = 0; i < seqRecords.size(); i++) {
                seqFile << i << "\t" << seqRecords[i].sequence << "\n";
            }
            seqFile.close();
        } else {
            std::cerr << "Error: Cannot open sequence file " << seqFileName << std::endl;
            continue;
        }

        int n = static_cast<int>(seqRecords.size());
        std::vector<std::vector<int>> distance_matrix(n, std::vector<int>(n, 0));

        Parameters params{1, -1, -2};

        #pragma omp parallel for collapse(2)
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                auto alignmentResult = needlemanWunsch(seqRecords[i].sequence, seqRecords[j].sequence, params);
                int d = std::max(0, 296 - alignmentResult.score);
                distance_matrix[i][j] = d;
                distance_matrix[j][i] = d;
            }
        }

        std::ofstream outputFile(outputFileName);
        if (outputFile.is_open()){
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    outputFile << distance_matrix[i][j] << " ";
                }
                outputFile << std::endl;
            }
            outputFile.close();

            runClustering(outputFileName, seqFileName);

        } else {
            std::cerr << "Error: Cannot open output file " << outputFileName << std::endl;
        }
    }

    return 0;
}