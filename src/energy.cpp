#include <armadillo>
#include <cassert>

#include "../include/energy.h"

double energy(arma::mat s) {
    assert (s.is_square() && "Spin configuration must be a square matrix");
    
    size_t L = s.n_rows;
    double E {};
    double E_ij {};

    int last = L-1;

    for (size_t i {}; i < L-1; i++) {
        for (size_t j {}; j < L-1; j++) {
            E_ij = s.at(i,j) * (s.at(i, j+1) + s.at(i+1, j));
            E -= E_ij;
        }

        size_t j {i};
        E_ij = s.at(i, last) * (s.at(i, 0) + s.at(i+1, last)) +
             s.at(last, j) * (s.at(0, j) + s.at(last, j+1));
        E -= E_ij;
    }
    
    E_ij = s.at(last, last) * (s.at(last, 0) + s.at(0, last));
    E -= E_ij;

    return E;
}