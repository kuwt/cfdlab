#include <iostream>
#include <string>

#include "Case.hpp"
#include <chrono>
int main(int argn, char **args) {

    if (argn > 1) {
        std::string file_name{args[1]};
        Case problem(file_name, argn, args);

        auto tstart = std::chrono::high_resolution_clock::now();
        problem.simulate();
        
        auto tstop = std::chrono::high_resolution_clock::now();
        auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(tstop - tstart);
        std::cout << "Time needed " << ms_int.count() << "ms\n";
    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
}
