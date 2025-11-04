#include <armadillo>

#include "../include/test.h"
#include "../include/energy.h"
#include "../include/heatcap_suscept.h"

bool test_energy1() {
    double tol {1e-12};

    arma::mat s {{-1, -1},
                 {-1, -1}};
    
    double E = energy(s);
    double E_expected = -8;
    double diff = std::abs(E - E_expected);

    return (diff < tol);
}

bool test_energy2() {
    double tol {1e-12};

    arma::mat s {{ 1, -1},
                 {-1, -1}};
    
    double E = energy(s);
    double E_expected = 0;
    double diff = std::abs(E - E_expected);

    return (diff < tol);    
}

bool test_energy3() {
    double tol {1e-12};

    arma::mat s {{ 1, 1},
                 {-1,-1}};
    
    double E = energy(s);
    double E_expected = 0;
    double diff = std::abs(E - E_expected);

    return (diff < tol);
}

bool test_energy4() {
    double tol {1e-12};

    arma::mat s {{ 1, -1},
                 { 1, -1}};

    double E = energy(s);
    double E_expected = 0;
    double diff = std::abs(E - E_expected);

    return (diff < tol);
}

bool test_energy5() {
    double tol {1e-12};

    arma::mat s {{ 1, -1},
                 {-1, 1}};

    double E = energy(s);
    double E_expected = 8;
    double diff = std::abs(E - E_expected);

    return (diff < tol);
}

bool test_energy6() {
    double tol {1e-12};

    arma::mat s {{1, -1},
                 {1,  1}};

    double E = energy(s);
    double E_expected = 0;
    double diff = std::abs(E - E_expected);

    return (diff < tol);
}

bool test_energy7() {
    double tol {1e-12};

    arma::mat s {{1, 1},
                 {1, 1}};

    double E = energy(s);
    double E_expected = -8;
    double diff = std::abs(E - E_expected);

    return (diff < tol);
}

bool test_energy8() {
    double tol {1e-12};

    arma::mat s {{ 1, -1, -1,  1},
                 { 1, -1,  1,  1},
                 {-1,  1, -1,  1},
                 { 1, -1, -1, -1}};
    
    double E = energy(s);
    double E_expected = 4;
    
    double diff = std::abs(E - E_expected);

    return (diff < tol);
}