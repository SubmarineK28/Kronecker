#include <iostream>
#include <iomanip>
#include <armadillo>
#include <complex>
#include <cmath>

using namespace arma;
using namespace std;

double exact_trace(const cx_mat& matrix, int power) {
    if (power == 0) {
        return matrix.n_rows;  
    }
    if (matrix.is_zero()) {
        return 0.0;
    }
    cx_mat matrix_power;
    try {
        matrix_power = arma::powmat(matrix, power);
    }
    catch (const std::runtime_error& e) {
        cerr << "Matrix power failed: " << e.what() << endl;
        return 0.0;
    }
    return trace(matrix_power).real();  
}
