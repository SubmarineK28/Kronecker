#include <iostream>
#include <iomanip>
#include <armadillo>
#include <complex>
#include <cmath>

using namespace arma;
using namespace std;

double Rayleigh_trace_estimator(const arma::cx_mat& A, int power_k, int M) {
    int n = A.n_rows;
    arma::vec trace_estimates(M, arma::fill::zeros);

    for (int i = 0; i < M; ++i) {
        arma::cx_vec z = arma::randn<arma::cx_vec>(n);
        z = arma::normalise(z) * std::sqrt(n);  
        arma::cx_vec Ak_z = z;
        for (int j = 0; j < power_k; ++j) {
            Ak_z = A * Ak_z;
        }
        trace_estimates(i) = arma::cdot(z, Ak_z).real();
    }

    return arma::mean(trace_estimates);
}
