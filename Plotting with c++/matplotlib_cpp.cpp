#include <vector>
#include <cmath>  // For sine function
#include "D:\C++ codes\matplotlibcpp\include\matplotlibcpp.h"

namespace plt = matplotlibcpp;

// Define M_PI if it's not already defined
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


// Run with:
// g++ "matplotlib_cpp.cpp" `
// >>  -I"C:\Users\mrsag\AppData\Local\Programs\Python\Python311\include" `
// >>  -I"C:\Users\mrsag\AppData\Local\Programs\Python\Python311\Lib\site-packages\numpy\core\include" `
// >>  -I"D:\C++ codes\matplotlibcpp\include" `
// >>  -L"C:\Users\mrsag\AppData\Local\Programs\Python\Python311\libs" `
// >>  -lpython311 -std=c++20 -mconsole -o matplotlib_cpp.exe
// >>

int main() {
    // Generate data for the sine wave
    std::vector<double> x, y;
    for (double i = -M_PI; i < 3 * M_PI; i += 0.1) {  // 2π range for one full sine wave
        x.push_back(i);                          // x values
        y.push_back(sin(i));                     // y values as sine of x
    }

    // Plot the sine wave
    plt::plot(x, y);
    plt::title("Sine Wave");
    plt::xlabel("x");
    plt::ylabel("sin(x)");
    plt::show();

    return 0;
}
