#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// Permet de contaténer une liste de plusieurs fichiers (fileNames) textes dans un nouveau fichier(outputFile)
void mergeFiles(const std::vector<std::string>& fileNames, const std::string& outputFile) {
    std::ofstream outputFileStream(outputFile); // nouveau fichier

    if (!outputFileStream.is_open()) {
        std::cerr << "Error: Could not open the output file: " << outputFile << '\n';
        return;
    }

    for (const auto& fileName : fileNames) {
        std::ifstream inputFile(fileName);

        if (!inputFile.is_open()) {
            std::cerr << "Error: Could not open input file: " << fileName << '\n';
            continue; // Skip to the next file if the current one cannot be opened
        }

        std::string line;
        while (std::getline(inputFile, line)) {
            outputFileStream << line << '\n'; // Write the line to the output file
        }

        inputFile.close();
    }

    // std::cout << "Files merged successfully into " << outputFile << "\n";
    outputFileStream.close();
}

