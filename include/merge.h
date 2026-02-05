#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// Allows you to contain a list of several text files (fileNames) in a new file (outputFile)
void mergeFiles(const std::vector<std::string>& fileNames, const std::string& outputFile) {
    std::ofstream outputFileStream(outputFile); // new file

    if (!outputFileStream.is_open()) {
        std::cerr << "Error: Could not open the output file: " << outputFile << '\n';
        return;
    }

    for (const auto& fileName : fileNames) {
        std::ifstream inputFile(fileName);

        if (!inputFile.is_open()) {
            std::cerr << "Error: Could not open input file: " << fileName << '\n';
            continue; // skip to the next file if the current one cannot be opened
        }

        std::string line;
        while (std::getline(inputFile, line)) {
            outputFileStream << line << '\n'; // write the line to the output file
        }

        inputFile.close();
    }

    outputFileStream.close();
}

