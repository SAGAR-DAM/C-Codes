#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;


// temp: 1

// void print_arr(std::vector<double> vect)
// {
//     for (int i = 0; i < vect.size(); ++i) 
//     {
//        cout<<vect[i]<<" ";
//     }
    
// }

// int main() {
//     int size;
//     std::cout << "Enter the size of the array: ";
//     std::cin >> size;

//     // Create a dynamic array using std::vector
//     std::vector<double> dynamicArray(size);

//     // Initialize array elements
//     for (int i = 0; i < size; ++i) {
//         dynamicArray[i] = pow(i,2);
//     }

//     print_arr(dynamicArray);
//     std::cout << std::endl;

//     return 0;
// }



// temp 2:


// int main() {
//     std::ofstream data("plot_data.dat");
//     for (float x = -5.0; x <= 5.0; x += 0.1) {
//         float y = sin(x);
//         data << x << " " << y << "\n";
//     }
//     data.close();

//     FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
//     fprintf(gnuplotPipe, "plot 'plot_data.dat' with lines\n");
//     fflush(gnuplotPipe);

//     std::cout << "Press enter to exit..." << std::endl;
//     std::cin.get();
//     return 0;
// }



// temp 3;
// #include <iostream>
// #include <fstream>
// #include <cmath>
// #include <cstdio> // For remove function

// int main() {
//     // Generate data for y = sin(x)
//     std::ofstream data1("data_sin.dat");
//     for (double x = -20.0; x <= 20.0; x += 0.01) {
//         double y = sin(x);
//         data1 << x << " " << y << "\n";
//     }
//     data1.close();

//     // Generate data for y = exp(cos(x))
//     std::ofstream data2("data_exp_cos.dat");
//     for (double x = -20.0; x <= 20.0; x += 0.01) {
//         double y = exp(cos(x));
//         data2 << x << " " << y << "\n";
//     }
//     data2.close();
    

//     // Open a pipe to GNUPLOT
//     FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
//     if (gnuplotPipe) {
//         // Send plotting commands to GNUPLOT
//         fprintf(gnuplotPipe, "set title 'Multiplot Example'\n");
//         fprintf(gnuplotPipe, "set xlabel 'x'\n");
//         fprintf(gnuplotPipe, "set ylabel 'y'\n");
//         fprintf(gnuplotPipe, "set grid\n");
//         fprintf(gnuplotPipe, "plot 'data_sin.dat' title 'y = sin(x)' with lines linecolor rgb 'green' linewidth 2, ");
//         fprintf(gnuplotPipe, "'data_exp_cos.dat' title 'y = exp(cos(x))' with linespoints pointtype 7 pointsize 1 linecolor rgb 'red' linewidth 2\n");
//         fflush(gnuplotPipe);

//         // Wait for the user to close the plot window
//         std::cout << "Press enter to exit..." << std::endl;
//         std::cin.get();

//         // Close the pipe to GNUPLOT
//         pclose(gnuplotPipe);

//         // Delete the data files
//         if (remove("data_sin.dat") != 0) {
//             std::cerr << "Error deleting data_sin.dat" << std::endl;
//         } else {
//             std::cout << "Deleted data_sin.dat" << std::endl;
//         }

//         if (remove("data_exp_cos.dat") != 0) {
//             std::cerr << "Error deleting data_exp_cos.dat" << std::endl;
//         } else {
//             std::cout << "Deleted data_exp_cos.dat" << std::endl;
//         }
//     } 
    
//     else {
//         std::cerr << "Error: Could not open pipe to GNUPLOT." << std::endl;
//     }

//     return 0;
// }

// temp 4:
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio> // For popen function

int main() {
    // Generate surface data
    std::ofstream data("surface_data.dat");
    for (double x = -5.0; x <= 5.0; x += 0.1) {
        for (double y = -5.0; y <= 5.0; y += 0.1) {
            double z = exp(-x*x-y*y)*sin(0.2*(x*x+y*y));
            data << x << " " << y << " " << z << "\n";
        }
        data << "\n";
    }
    data.close();

    // Write GNUPLOT script
    std::ofstream gp_script("plot_surface.gp");
    gp_script << "set title 'Interactive Surface Plot'\n";
    gp_script << "set xlabel 'x'\n";
    gp_script << "set ylabel 'y'\n";
    gp_script << "set zlabel 'z'\n";
    gp_script << "set pm3d\n";
    gp_script << "set hidden3d\n";
    gp_script << "splot 'surface_data.dat' with lines\n";
    gp_script.close();

    // Execute GNUPLOT script
    FILE *gnuplotPipe = popen("gnuplot", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "load 'plot_surface.gp'\n");
        fflush(gnuplotPipe);

        std::cout << "Press enter to exit..." << std::endl;
        std::cin.get();

        pclose(gnuplotPipe);
    } else {
        std::cerr << "Error: Could not open pipe to GNUPLOT." << std::endl;
    }

    // Optionally, delete the data file
    if (remove("surface_data.dat") != 0) {
        std::cerr << "Error deleting surface_data.dat" << std::endl;
    } else {
        std::cout << "Deleted surface_data.dat" << std::endl;
    }

    return 0;
}
