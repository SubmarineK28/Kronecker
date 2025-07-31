#include <iostream>
#include <iomanip>
#include <armadillo>
#include <complex>
#include <cmath>

using namespace arma;
using namespace std;

double Hutchinson_trace_estimator(const arma::cx_mat& A, int power_k, int M) {
    int n = A.n_rows;
    arma::vec trace_estimates(M, arma::fill::zeros);

    for (int i = 0; i < M; ++i) {
        arma::cx_vec z = arma::randi<arma::cx_vec>(n, arma::distr_param(0, 1));
        z = 2 * z - 1; 
        
        arma::cx_vec Ak_z = z;
        for (int j = 0; j < power_k; ++j) {
            Ak_z = A * Ak_z;
        }

        trace_estimates(i) = arma::cdot(z, Ak_z).real();
    }

    return arma::mean(trace_estimates);
}
