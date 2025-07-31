
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <complex>
#include <cmath>

using namespace arma;
using namespace std;

cx_mat createPauliX() {
    cx_mat X(2, 2);
    X(0, 0) = cx_double(0, 0);
    X(0, 1) = cx_double(1, 0);
    X(1, 0) = cx_double(1, 0);
    X(1, 1) = cx_double(0, 0);
    return X;
}

cx_mat createPauliY() {
    cx_mat Y(2, 2);
    Y(0, 0) = cx_double(0, 0);
    Y(0, 1) = cx_double(0, -1);
    Y(1, 0) = cx_double(0, 1);
    Y(1, 1) = cx_double(0, 0);
    return Y;
}

cx_mat createPauliZ() {
    cx_mat Z(2, 2);
    Z(0, 0) = cx_double(1, 0);
    Z(0, 1) = cx_double(0, 0);
    Z(1, 0) = cx_double(0, 0);
    Z(1, 1) = cx_double(-1, 0);
    return Z;
}

arma::cx_vec apply_H0(const arma::cx_vec& z, int N) {
    int dim = 1 << N;  // 2^N
    arma::cx_vec result(dim, arma::fill::zeros);

    for (int j = 0; j < N; ++j) {
        for (int k = j + 1; k < N; ++k) {
            double factor = 1.0 / std::pow(k - j, 3);

            int bit_j = N - 1 - j;
            int bit_k = N - 1 - k;

            for (int state = 0; state < dim; ++state) {
                int bj = (state >> bit_j) & 1;
                int bk = (state >> bit_k) & 1;

                // --------- ?^x_j ?^x_k and ?^y_j ?^y_k --------------
                int flipped = state ^ (1 << bit_j) ^ (1 << bit_k);

                // ?^x ? ?^x: always contributes 1
                result[state] += factor * z[flipped];

                // ?^y ? ?^y: contributes +1 or -1 depending on bj, bk
                std::complex<double> phase = 0.0;
                if (bj == bk) phase = -1.0;
                else phase = 1.0;

                result[state] += factor * phase * z[flipped];

                // --------- -2 * ?^z_j ?^z_k --------------
                double zj = (bj == 0) ? 1.0 : -1.0;
                double zk = (bk == 0) ? 1.0 : -1.0;
                result[state] += (-2.0 * factor) * zj * zk * z[state];
            }
        }
    }

    return result;
}

