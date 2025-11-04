#pragma once 
#include <armadillo>

/**
 * @brief Performs a Monte Carlo cycle on a given spin configuration
 * @param s Spin configuration as a square matrix
 * @param L Number of rows/columns in lattice
 * @param E_temp Current energy of the spin configuration
 * @param Esq_temp Current energy squared of the spin configuration
 * @param M_temp Current magnetisation of the spin configuration
 * @param absM_temp Current absolute value of magnetisation of the spin configuration
 * @param Msq_temp Current magnetisation squared of the spin configuration
 * @param E_sum Current sum of energies of the spin configurations encountered during Monte Carlo steps so far
 * @param Esq_sum Current sum of energies squared of the spin configurations encountered during Monte Carlo steps so far
 * @param absM_sum Current sum of absolute value of magnetisation of the spin configurations encountered during Monte Carlo steps so far
 * @param Msq_sum Current sum of magnetisation squared of the spin configurations encountered during Monte Carlo steps so far
 * @param w4 Probability of accepting a spin flip when the energy difference is 4
 * @param w8 Probability of accepting a spin flip when the energy difference is 8
 * @param generator Random number generator
 * @param filename_ϵ Name of the file to write the energy per spin to. If empty, no file is written
 * @param filename_m Name of the file to write the magnetisation per spin to. If empty, no file is written
 * @return none
*/
void monte_carlo_cycle(arma::mat &s, double &E_temp, double &Esq_temp, double &M_temp, double &absM_temp, double &Msq_temp, double &E_sum, double &Esq_sum, double &absM_sum, double &Msq_sum, double w4, double w8, std::mt19937 &generator, std::string filename_ϵ = "", std::string filename_m = "");
