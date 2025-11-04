#pragma once
#include <armadillo>
/**
 * @brief Calculates the energy of a given spin configuration.
 * @param s Spin configuration as a square matrix
 * @param L Number of rows/columns in lattice
 * @return The total energy E of the spin configuration 
*/
double energy(arma::mat s);
