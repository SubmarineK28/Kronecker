#include <iostream>
#include <iomanip>
#include <armadillo>
#include <complex>
#include <cmath>

using namespace arma;
using namespace std;

double Unit_vector_estimator(const arma::cx_mat& H, int k, int M) {
    int n = H.n_rows;
    arma::vec trace_estimates(M, arma::fill::zeros);

    for (int i = 0; i < M; ++i) {
        arma::cx_vec z = arma::zeros<arma::cx_vec>(n);
        int index = arma::randi<int>(arma::distr_param(0, n - 1));
        z(index) = arma::cx_double(1,0);

        arma::cx_vec Hk_z = z;
        for (int j = 0; j < k; ++j) {
            Hk_z = H * Hk_z; 
        }
        trace_estimates(i) = arma::cdot(z, Hk_z).real();
    }
    return arma::mean(trace_estimates);
}
