#include "check_1d.hpp"
#include <exception>
#include <iostream>
int main() {
    try {
        euler::check::all();
        std::cout << "All 1D checks passed\n";
    } catch (const std::exception &error) {
        std::cerr << "CHECK FAILED: " << error.what() << '\n';
        return 1;
    }
}
