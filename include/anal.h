#pragma once
#include <tuple>

/**
 * @brief Calculates the analytical expectation values for the energy of a 2x2 lattice
 * @param β Inverse temperature
 * @return Tuple containing the analytical expectation values for the energy per spin ϵ and the energy per spin squared ϵ²
*/
std::tuple<double,double> ϵ_anal(double β);

/**
 * @brief Calculates the analytical expectation values for the magnetisation per spiun of a 2x2 lattice
 * @param β Inverse temperature
 * @return Tuple containing the analytical expectation values for the absolute value of the magnetisation per spin m and the magnetisation per spin squared m²
*/
std::tuple<double,double> m_anal(double β);

/**
 * @brief Calculates the analytical value for the heat capacity normalised to the number of spins of a 2x2 lattice
 * @param β Inverse temperature
 * @return The analytical value for the heat capacity normalised to the number of spins C_V
*/
double C_V_anal(double β);

/**
 * @brief Calculates the analytical value for the susceptibility normalised to the number of spins of a 2x2 lattice
 * @param β Inverse temperature
 * @return The analytical value for the susceptibility normalised to the number of spins χ
*/
double χ_anal(double β);