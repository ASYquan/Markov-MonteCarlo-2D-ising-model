#include <iostream>
#include <random>
#include <cmath>
#include <cassert>
#include <chrono>
#include <armadillo>
#include "../include/MCMC.h"

void monte_carlo_cycle(arma::mat &s, double &E_temp, double &Esq_temp, double &M_temp, double &absM_temp, double &Msq_temp, double &E_sum, double &Esq_sum, double &absM_sum, double &Msq_sum, double w4, double w8, std::mt19937 &generator, std::string filename_ϵ, std::string filename_m) {
    assert(s.is_square() && "Spin configuration must be a square matrix");
    assert(w4 >= 0.0 && w4 <= 1.0 && "Probability of accepting a spin flip when the energy difference is 4 must be between 0 and 1");
    assert(w8 >= 0.0 && w8 <= 1.0 && "Probability of accepting a spin flip when the energy difference is 8 must be between 0 and 1");

    // Random number generator
    // Uniform distribution
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    size_t L = s.n_rows;
    int n_spins = L * L;
    double n_spins_d = static_cast<double>(n_spins);

    arma::mat ϵ(L, L);
    arma::mat m(L, L);

    for (size_t i = 0; i < L; ++i) {
        for (size_t j = 0; j < L; ++j) {
            // Calculate the energy difference
            int up    = (i == 0)     ? (L - 1) : (i - 1);
            int down  = (i == L - 1) ? 0       : (i + 1);
            int left  = (j == 0)     ? (L - 1) : (j - 1);
            int right = (j == L - 1) ? 0       : (j + 1);

            double ΔE = 2.0 * s.at(i, j) * (s.at(up, j) + s.at(down, j) + s.at(i, left) + s.at(i, right));
            double ΔM = - 2.0 * s.at(i, j);

            // Flip the spin if the energy difference is negative or zero
            if (ΔE <= 0) {
                s.at(i, j) *= -1; // Flip the spin

                // Update current energy and magnetisation values
                E_temp += ΔE;
                Esq_temp = E_temp * E_temp;
                M_temp += ΔM;
                absM_temp = std::abs(M_temp);
                Msq_temp = M_temp * M_temp;
            } 
            
            else {
                double w = (ΔE == 4) ? w4 : ((ΔE == 8) ? w8 : 0.0);

                // Flip the spin if w is larger than or equal to a random number r
                if (dis(generator) <= w) {
                    s.at(i, j) *= -1; // Flip the spin

                    // Update current energy and magnetisation values
                    E_temp += ΔE;
                    Esq_temp = E_temp * E_temp;
                    M_temp += ΔM;
                    absM_temp = std::abs(M_temp);
                    Msq_temp = M_temp * M_temp;
                }
            }

            // Update sums over energies and magnetisation values
            E_sum += E_temp;
            Esq_sum += Esq_temp;
            absM_sum += absM_temp;
            Msq_sum += Msq_temp;

            // Storing the energy and magnetisation
            ϵ.at(i, j) = E_temp / n_spins_d;
            m.at(i, j) = M_temp / n_spins_d;
        }
    }

    // Add the stored energy values to file
    if (!filename_ϵ.empty()) {
        std::ofstream file(filename_ϵ, std::ios::app);
        arma::vectorise(ϵ).raw_print(file);
    }

    // Add the stored magnetisation values to file
    if (!filename_m.empty()) {
        std::ofstream file(filename_m, std::ios::app);
        arma::vectorise(m).raw_print(file);
    }
}
