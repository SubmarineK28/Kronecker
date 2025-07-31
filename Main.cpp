#include <iostream>
#include <iomanip>
#include <armadillo>
#include <complex>
#include <cmath>

using namespace arma;
using namespace std;

cx_mat createH0(int N);
cx_mat createM0(int N);
cx_mat createFullHamiltonian(int N, double beta);

//double exact_trace(const cx_mat& matrix, int power);
//double Gaussian_trace_estimator(const arma::cx_mat& H, int k, int M);
//double Rayleigh_trace_estimator(const arma::cx_mat& A, int power_k, int num_samples);
//double Hutchinson_trace_estimator(const arma::cx_mat& A, int power_k, int M);
//double Unit_vector_estimator(const arma::cx_mat& H, int k, int M);

double Rayleigh_trace_estimator_matrix_free(int N, int power_k, int M);

//int main() {
//    int k = 1, M = 1000;
//    int N = 2;
//    double beta = 0;
//
//    cout << setw(21) << "Gauss_tr"
//        << setw(12) << "Ray_tr"
//        << setw(12) << "Hutch_tr"
//        << setw(12) << "Un_vec" << "\n";

    //for (int i = 0; i < 8; ++i) {

        //cx_mat H0 = createH0(N);
        //cout << "H0:\n" << H0 << endl;
        //cx_mat M0 = createM0(N);
        //cout << "M0:\n" << M0 << endl;
        //cx_mat H = createFullHamiltonian(N, beta);
        //cout << "Exact Trace of H = " << exact_tr << "\n";

        //double exact_tr = exact_trace(H, k);
        //double Gauss_tr = Gaussian_trace_estimator(H, k, M);
        //double Ray_tr = Rayleigh_trace_estimator(H, k, M);
        //double Hutch_tr = Hutchinson_trace_estimator(H, k, M);
        //double Un_vec = Unit_vector_estimator(H, k, M);

        //cout << fixed << setprecision(6);
        //cout << "N = " << setw(2) << N << " | "
        //    << setw(12) << Gauss_tr
        //    << setw(12) << Ray_tr
        //    << setw(12) << Hutch_tr
        //    << setw(12) << Un_vec
        //    << endl;

        //N += 1;
    //}

    //double trace_approx = Rayleigh_trace_estimator_matrix_free(N, k, M);
    //std::cout << "Trace ? " << trace_approx << std::endl;

    //return 0;
//}


int main() {
    int k = 1, M = 1000;
    int N = 2;
    double beta = 0;

    cout << setw(21) << "Gauss_tr"
        << setw(12) << "Ray_tr"
        << setw(12) << "Hutch_tr"
        << setw(12) << "Un_vec" << "\n";

    double trace_approx = Rayleigh_trace_estimator_matrix_free(N, k, M);
    std::cout << "Trace ? " << trace_approx << std::endl;

    return 0;
}