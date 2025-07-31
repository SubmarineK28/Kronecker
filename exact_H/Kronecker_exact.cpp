#include <iostream>
#include <iomanip>
#include <armadillo>
#include <complex>
#include <cmath>

using namespace arma;
using namespace std;

// --- Паули матрицы с комплексными числами ---
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

cx_mat createH0(int N) {
    cx_mat H0;
    bool isFirstTerm = true;

    for (int j = 1; j <= N; ++j) {
        for (int k = j + 1; k < N; ++k) {
            double factor = 1.0 / pow(k - j, 3);

            for (const auto& sigma_c : { createPauliX(), createPauliY() }) {
                cx_mat term = eye<cx_mat>(1, 1);

                for (int pos = 0; pos < N; ++pos) {
                    if (pos == j || pos == k) {
                        term = kron(term, sigma_c);
                    }
                    else {
                        term = kron(term, eye<cx_mat>(2, 2));
                    }
                }

                if (isFirstTerm) {
                    H0 = factor * term;
                    isFirstTerm = false;
                }
                else {
                    H0 += factor * term;
                }
            }

            // σz ⊗ σz
            cx_mat termZ = eye<cx_mat>(1, 1);
            for (int pos = 0; pos < N; ++pos) {
                if (pos == j || pos == k) {
                    termZ = kron(termZ, createPauliZ());
                }
                else {
                    termZ = kron(termZ, eye<cx_mat>(2, 2));
                }
            }

            H0 += (-2 * factor) * termZ;
        }
    }

    return H0;
}

cx_mat createM0(int N) {
    cx_mat M0 = zeros<cx_mat>(pow(2, N), pow(2, N));

    for (int j = 0; j < N; ++j) {
        cx_mat term = eye<cx_mat>(1, 1);  // Начинаем с 1x1 единичной матрицы

        for (int pos = 0; pos < N; ++pos) {
            if (pos == j) {
                term = kron(term, createPauliZ());
            }
            else {
                term = kron(term, eye<cx_mat>(2, 2));
            }
        }

        M0 += term;
    }

    return M0;
}

cx_mat createFullHamiltonian(int N, double beta) {
    cx_mat H0 = createH0(N);
    cx_mat M0 = createM0(N);
    return H0 - beta * M0;
}



