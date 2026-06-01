#include "clustering.h"
#include <iostream>
#include <cstdlib>

// function to run the clustering Python script with the generated distance matrix and sequences files
void runClustering(const std::string& distanceMatrixFile, const std::string& sequencesFile) {
    std::string command = "python3 clustering.py --distance_matrix " + distanceMatrixFile + " --sequences " + sequencesFile;
    int result = system(command.c_str());
    if (result != 0) {
        std::cerr << "Error executing clustering script." << std::endl;
    }
}
