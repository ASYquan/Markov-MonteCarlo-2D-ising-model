#pragma once

/**
 * @brief Calculates the heat capacity normalised to the number of spins
 * @param n_spins Number of spins
 * @param β Inverse temperature
 * @param exp_ϵ Expectation value for the energy per spin
 * @param exp_ϵsq Expectation value for the energy per spin squared
 * @return The heat capacity normalised to the number of spins C_V
*/
double C_V(std::size_t n_spins, double β, double exp_ϵ, double exp_ϵsq);

/**
 * @brief Calculates the susceptibility normalised to the number of spins
 * @param n_spins Number of spins
 * @param β Inverse temperature
 * @param exp_absm Expectation value for the magnetisation per spin
 * @param exp_msq Expectation value for the magnetisation per spin squared
 * @return The susceptibility normalised to the number of spins χ
*/
double χ(std::size_t n_spins, double β, double exp_absm, double exp_msq);
