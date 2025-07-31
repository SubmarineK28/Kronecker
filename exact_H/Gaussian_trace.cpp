
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <complex>
#include <cmath>

using namespace arma;
using namespace std;

double Gaussian_trace_estimator(const arma::cx_mat& H, int k, int M) {
    int n = H.n_rows;
    arma::vec trace_estimates(M, arma::fill::zeros);

    for (int i = 0; i < M; ++i) {
        arma::cx_vec z = arma::randn<arma::cx_vec>(n);

        arma::cx_vec Hk_z = z;
        for (int j = 0; j < k; ++j) {
            Hk_z = H * Hk_z; 

        trace_estimates(i) = arma::cdot(z, Hk_z).real();
    }

    return arma::mean(trace_estimates);
}
