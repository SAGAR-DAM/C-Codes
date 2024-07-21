#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

// Function to check if a number is prime
bool is_prime(int num) {
    if (num <= 1) return false;
    if (num <= 3) return true;
    if (num % 2 == 0 || num % 3 == 0) return false;

    for (int i = 5; i * i <= num; i += 6) {
        if (num % i == 0 || num % (i + 2) == 0) return false;
    }

    return true;
}

// Function to find the nth prime
int nth_prime(int n) {
    if (n < 1) {
        std::cerr << "n must be a positive integer." << std::endl;
        return -1;
    }

    int count = 0;
    int num = 1;

    while (count < n) {
        num++;
        if (is_prime(num)) {
            count++;
        }
    }

    return num;
}

int main() {
    int n;
    std::cout << "Enter the value of n: ";
    std::cin >> n;
    std::cout << "\n\n";

    std::vector<int> primes;
    std::vector<int> indices;

    for (int j = 1; j <= n; j++) {
        int result = nth_prime(j);
        if (result != -1) {
            indices.push_back(j);
            primes.push_back(result);
        }
    }

    // Write data to file
    std::ofstream data_file("data_prime.dat");
    for (size_t i = 0; i < primes.size(); i++) {
        data_file << indices[i] << " " << primes[i] << "\n";
    }
    data_file.close();

    // Open a pipe to GNUPLOT
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        // Send plotting commands to GNUPLOT
        fprintf(gnuplotPipe, "set title 'Index vs Prime Number (up to the %dth prime)'\n", n);
        fprintf(gnuplotPipe, "set xlabel 'Index (i)'\n");
        fprintf(gnuplotPipe, "set ylabel 'Prime Number'\n");
        fprintf(gnuplotPipe, "set grid\n");
        fprintf(gnuplotPipe, "plot 'data_prime.dat' with linespoints pointtype 7 pointsize 1 linecolor rgb 'red' title 'Prime Numbers'\n");
        fflush(gnuplotPipe);

        // Wait for the user to close the plot window
        std::cout << "Press enter to delete the data file and exit..." << std::endl;
        std::cin.get();  // This will consume the newline character left in the input buffer
        std::cin.get();  // This will wait for the actual Enter key press

        // Close the pipe to GNUPLOT
        pclose(gnuplotPipe);

        // Delete the data file
        if (remove("data_prime.dat") != 0) {
            std::cerr << "Error deleting data_prime.dat" << std::endl;
        } else {
            std::cout << "Deleted data_prime.dat" << std::endl;
        }
    } else {
        std::cerr << "Error: Could not open pipe to GNUPLOT." << std::endl;
    }

    return 0;
}
