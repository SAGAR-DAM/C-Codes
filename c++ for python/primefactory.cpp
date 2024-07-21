#include <iostream>
#include <cmath>
#include <vector>

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
    std::cin >> n;

    std::vector<int> primes;
    for (int j = 1; j <= n; j++) {
        int result = nth_prime(j);
        if (result != -1) {
            primes.push_back(result);
        }
    }

    // Output the primes
    for (int i = 0; i < primes.size(); ++i) {
        std::cout << primes[i];
        if (i != primes.size() - 1) {
            std::cout << " ";
        }
    }

    return 0;
}
