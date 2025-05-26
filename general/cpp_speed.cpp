#include <iostream>
#include <chrono>

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    long long x = 100;
    for (long long i = 0; i < 1000000000LL; ++i) 
    {
        x=x+1;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout<<"x: "<<x<<std::endl;
    std::cout << "Time taken: " << elapsed.count() << " seconds\n";
    return 0;
}
