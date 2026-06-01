#include "parse.h"
#include "algorithm.h"
#include "clustering.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <numeric>

constexpr int MINIMUM_CLUSTER_SIZE = 20;
constexpr int DISTANCE_THRESHOLD = 10;


int main() {
    namespace fs = std::filesystem;
    // path to the directory containing the Fastq files
    const std::string dirPath = "Bioinformatika - jeleni-2/fastq";

    const fs::path matricesDir = "matrices";
    const fs::path sequencesDir = "sequences";
    const fs::path clusteringDir = "clustering";

    // create output directories if they don't exist
    fs::create_directories(matricesDir);
    fs::create_directories(sequencesDir);

    // loop through all Fastq files in the specified directory
    
    for (const auto& entry : fs::directory_iterator(dirPath)) {
        if (!entry.is_regular_file()) continue;
        if (entry.path().extension() != ".fastq") continue;
        std::string baseName = entry.path().filename().string();
        if (baseName.empty() || baseName[0] != 'J') continue;

        std::string filename = entry.path().string();
        std::cout << "Processing file: " << entry.path().filename().string() << std::endl;

        auto result = parseFastq(filename);
        auto allRecords = result.first;
        auto seqRecords = result.second;

        if (seqRecords.empty()) {
            std::cerr << "Warning: no sequences parsed from " << filename << std::endl;
            continue;
        }

        // derive sample name from the file name
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
        // initialize distance matrix with zeros
        std::vector<std::vector<int>> distance_matrix(n, std::vector<int>(n, 0));

        // alignment parameters: match score = 1, mismatch penalty = -1, gap penalty = -2
        Parameters params{1, -1, -2};

        // compute pairwise needleman-wunsch distances in parallel using OpenMP
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                auto alignmentResult = needlemanWunsch(seqRecords[i].sequence, seqRecords[j].sequence, params);
                int d = std::max(0, 296 - alignmentResult.score);
                distance_matrix[i][j] = d;
                distance_matrix[j][i] = d;
            }
        }

        // write distance matrix to output file
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
    

    std::vector<std::string> representatives;

    // loop through all clustering output files to extract representatives of clusters larger than the minimum size
    for (const auto& entry : fs::directory_iterator(clusteringDir)) {
        if (!entry.is_regular_file()) continue;

        std::ifstream clusteringFile(entry.path());
        if (!clusteringFile.is_open()) {
            std::cerr << "Error: Cannot open clustering file " << entry.path().string() << std::endl;
            continue;
        }

        std::string line;
        while (std::getline(clusteringFile, line)) {
            const std::string sizeTag = "size: ";
            const std::string consensusTag = "consensus: ";

            auto sizePos = line.find(sizeTag);
            auto consensusPos = line.find(consensusTag);
            if (sizePos == std::string::npos || consensusPos == std::string::npos) continue;

            int clusterSize = std::stoi(line.substr(sizePos + sizeTag.size(), consensusPos - (sizePos + sizeTag.size())));
            if (clusterSize <= MINIMUM_CLUSTER_SIZE) continue;

            representatives.push_back(line.substr(consensusPos + consensusTag.size()));
        }
    }

    if (representatives.empty()) {
        std::cout << "No representatives found for grouping." << std::endl;
        return 0;
    }

    const Parameters params{1, -1, -2};
    const int distanceThreshold = DISTANCE_THRESHOLD;
    const int n = static_cast<int>(representatives.size());

    // compute pairwise distances between representatives using the same distance function as before
    auto sequenceDistance = [&](const std::string& left, const std::string& right) -> int {
        auto alignmentResult = needlemanWunsch(left, right, params);
        return std::max(0, 296 - alignmentResult.score);
    };

    std::vector<int> parent(n);
    std::iota(parent.begin(), parent.end(), 0);

    // perform union-find to group representatives based on distance threshold
    auto findRoot = [&](auto&& self, int node) -> int {
        if (parent[node] == node) return node;
        parent[node] = self(self, parent[node]);
        return parent[node];
    };

    auto unite = [&](int left, int right) {
        int rootLeft = findRoot(findRoot, left);
        int rootRight = findRoot(findRoot, right);
        if (rootLeft != rootRight) {
            parent[rootRight] = rootLeft;
        }
    };

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (sequenceDistance(representatives[i], representatives[j]) < distanceThreshold) {
                unite(i, j);
            }
        }
    }

    std::unordered_map<int, std::vector<std::string>> groupedRepresentatives;
    for (int i = 0; i < n; ++i) {
        groupedRepresentatives[findRoot(findRoot, i)].push_back(representatives[i]);
    }

    fs::path groupedOutputPath = clusteringDir / "representative_groups.txt";
    std::ofstream groupedOutput(groupedOutputPath);
    if (!groupedOutput.is_open()) {
        std::cerr << "Error: Cannot open grouped output file " << groupedOutputPath.string() << std::endl;
        return 0;
    }

    // for each group of representatives, find the one with the smallest total distance to all other members and write the group info to the output file
    int groupIndex = 0;
    for (const auto& [root, group] : groupedRepresentatives) {
        std::string groupRepresentative = group.front();
        int bestTotalDistance = std::numeric_limits<int>::max();

        for (const auto& candidate : group) {
            int totalDistance = 0;
            for (const auto& other : group) {
                if (candidate == other) continue;
                totalDistance += sequenceDistance(candidate, other);
            }

            if (totalDistance < bestTotalDistance) {
                bestTotalDistance = totalDistance;
                groupRepresentative = candidate;
            }
        }

        groupedOutput << "Group " << groupIndex++
                      << " | size: " << group.size()
                      << " | representative: " << groupRepresentative << "\n";
        for (const auto& member : group) {
            groupedOutput << member << "\n";
        }
        groupedOutput << "\n";
    }

    std::cout << "Stored " << groupedRepresentatives.size()
              << " representative groups in " << groupedOutputPath.string() << std::endl;


    return 0;
}