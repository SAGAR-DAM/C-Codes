#include <iostream>
#include <iomanip> // For std::setw
#include <string>  // For std::string

// ANSI color codes
const std::string GREEN = "\033[92m";
const std::string RESET = "\033[0m";

// Function to display the progress bar
void display_progress_bar(int total_iterations) {
    std::cout << "Processing" << std::string(100 - 10, ' ') << std::flush; // Title + space for bar

    for (int i = 0; i < total_iterations; ++i) {
        // Calculate percentage and number of hashes for the bar
        float percentage = (i + 1) / static_cast<float>(total_iterations) * 100;
        int bar_width = 70; // Width of the progress bar in characters
        int filled_width = static_cast<int>(bar_width * percentage / 100);
        
        // Construct the bar
        std::string bar = GREEN + std::string(filled_width, '#') + std::string(bar_width - filled_width, '-') + RESET;

        // Output progress bar
        std::cout << "\r" << bar << "| " << std::fixed << std::setprecision(0) << percentage 
                  << "% " << std::flush; // '\r' to return to the start of the line
    }

    // Finish the bar
    std::cout << std::endl;
}

int main() {
    const int total_iterations = 100000; // You can change this to any value
    display_progress_bar(total_iterations);
    return 0;
}
