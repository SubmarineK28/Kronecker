#include <iostream>
#include <iomanip>
#include <armadillo>
#include <complex>
#include <cmath>

using namespace arma;
using namespace std;

arma::cx_vec apply_H0(const arma::cx_vec& z, int N);

double Rayleigh_trace_estimator_matrix_free(int N, int power_k, int M) {
    int dim = 1 << N;
    arma::vec trace_estimates(M, arma::fill::zeros);

    for (int i = 0; i < M; ++i) {
        arma::cx_vec z = arma::randn<arma::cx_vec>(dim);
        z = arma::normalise(z) * std::sqrt(dim);

        arma::cx_vec Ak_z = z;
        for (int j = 0; j < power_k; ++j) {
            Ak_z = apply_H0(Ak_z, N);
        }

        trace_estimates(i) = arma::cdot(z, Ak_z).real();
    }

    return arma::mean(trace_estimates);
}
