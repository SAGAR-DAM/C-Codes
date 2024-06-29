#include <iostream>
#include <random>
#include <vector>
#include <iomanip> // For std::setprecision

int main() {
    // Create a random device and a generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Define the distribution to generate floating point numbers between 0 and 1
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Define the size of the array
    const int array_size = 50;
    std::vector<double> random_array(array_size);

    // Fill the array with random floating point numbers between 0 and 1
    for (int i = 0; i < array_size; ++i) {
        random_array[i] = dis(gen);
    }

    // Print the array
    std::cout << "Random array of size " << array_size << " with floating point numbers between 0 and 1:" << std::endl;
    std::cout << std::fixed << std::setprecision(4); // Set precision for printing
    for (int i = 0; i < array_size; ++i) {
        std::cout << random_array[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
