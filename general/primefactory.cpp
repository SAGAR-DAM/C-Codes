#include <iostream>
#include <cmath>

// Function to check if a number is prime
bool isPrime(int num) {
    if (num <= 1) return false;
    if (num <= 3) return true;
    if (num % 2 == 0 || num % 3 == 0) return false;

    for (int i = 5; i * i <= num; i += 6) {
        if (num % i == 0 || num % (i + 2) == 0) return false;
    }

    return true;
}

// Function to find the nth prime
int nthPrime(int n) {
    if (n < 1) {
        std::cerr << "n must be a positive integer." << std::endl;
        return -1;
    }

    int count = 0;
    int num = 1;

    while (count < n) {
        num++;
        if (isPrime(num)) {
            count++;
        }
    }

    return num;
}

int main() {
    int n;
    std::cout << "Enter the value of n: ";
    std::cin >> n;
    std::cout<<"\n\n";
    int result;

    for(int j = 1; j<=n ; j++)
    {
        result = nthPrime(j);
        if (result != -1) 
        {
            std::cout << "The " << j << "th prime number is " << result << std::endl;
        }
    }

    return 0;
}
