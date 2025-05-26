#include <iostream>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

// Function to recursively print directory structure like the 'tree' command
void printDirectoryTree(const fs::path& dirPath, const std::string& prefix = "") {
    try {
        if (!fs::exists(dirPath)) {
            std::cerr << "Directory path " << dirPath << " does not exist.\n";
            return;
        }
        if (!fs::is_directory(dirPath)) {
            std::cerr << dirPath << " is not a directory.\n";
            return;
        }

        std::cout << "Processing directory: " << dirPath << std::endl; // Debug output

        // Iterate over the directory contents
        int count = 0;
        for (const auto& entry : fs::directory_iterator(dirPath)) {
            count++;
            bool isLast = (std::next(fs::directory_iterator(entry)) == fs::directory_iterator());

            std::cout << prefix;
            if (isLast) {
                std::cout << "\\-- ";  // Use ASCII \-- for last entry
            } else {
                std::cout << "|-- ";   // Use ASCII |-- for non-last entries
            }

            std::cout << entry.path().filename().string() << "\n";

            // If the entry is a directory, recursively print its contents
            if (fs::is_directory(entry.status())) {
                std::string newPrefix = prefix + (isLast ? "    " : "|   ");
                printDirectoryTree(entry.path(), newPrefix);
            }
        }

        if (count == 0) {
            std::cerr << "No files or directories found in " << dirPath << "\n";
        }

    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error accessing directory: " << e.what() << "\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <directory-path>\n";
        return 1;
    }

    fs::path inputPath = argv[1];

    std::cout << "Received input path: " << inputPath << std::endl; // Debug output

    // Call the function to print the directory structure
    printDirectoryTree(inputPath);

    return 0;
}
