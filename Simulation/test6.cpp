#include <iostream>
#include <vector>
#include <initializer_list>

// Concatenate function
std::vector<double> concatenate(const std::initializer_list<std::vector<double>> &vectors)
{
    std::vector<double> result;

    // Reserve total size to avoid reallocations
    size_t total_size = 0;
    for (const auto &vec : vectors)
        total_size += vec.size();
    result.reserve(total_size);

    // Insert all elements
    for (const auto &vec : vectors)
        result.insert(result.end(), vec.begin(), vec.end());

    return result;
}

int main()
{
    // Define input vectors
    std::vector<double> x = {1.0, 2.0};
    std::vector<double> y = {3.0, 4.0, 5.0};
    std::vector<double> z = {6.0};

    // Merge vectors
    std::vector<double> merged = concatenate({x, y, z});

    // Print result
    std::cout << "Merged vector: ";
    for (const auto &val : merged)
        std::cout << val << " ";
    std::cout << std::endl;

    return 0;
}
