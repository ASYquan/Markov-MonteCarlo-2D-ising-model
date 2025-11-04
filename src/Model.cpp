#include <iostream>
#include <armadillo>
#include <random>
#include <tuple>
#include <chrono>
#include <stdio.h>  
#include <omp.h>
#include "../include/energy.h"
#include "../include/heatcap_suscept.h"
#include "../include/MCMC.h"
#include "../include/Model.h"

Model::Model(std::size_t L, double T, std::size_t n_checks, std::size_t n_cycles_per_check, int seed) 
            : L(L), n_checks(n_checks), n_cycles_per_check(n_cycles_per_check), seed(seed) {}

Model::~Model() {}

// void Model::MCMC(bool parallel = false, bool ordered = false, std::string filename = "") {


//     arma::arma_rng::set_seed(seed);
//     std::mt19937 generator(seed);

//     if (ordered) {
//         arma::mat s = arma::mat(L, L, arma::fill::ones);
//     }
//     else {
//         arma::mat s = arma::mat(L, L, arma::fill::randu);
//         s = arma::sign(s - 0.5);
//     }

//     double β = 1 / T;
//     double w4 = std::exp(- β * 4);
//     double w8 = std::exp(- β * 8);
//     double E_temp = energy(s);   
//     double Esq_temp = E_temp * E_temp;  
//     double M_temp = arma::accu(s);
//     double absM_temp = std::abs(M_temp);   
//     double Msq_temp = M_temp * M_temp;     
//     double E_sum = E_temp;    
//     double Esq_sum = Esq_temp;  
//     double absM_sum = absM_temp;   
//     double Msq_sum = Msq_temp;  

//     for (size_t i = 0; i < n_checks; ++i) {
//         for (size_t j = 0; j < n_cycles_per_check; ++j) {
//             monte_carlo_cycle(s, generator, w4, w8, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp);
//         }
//         double i_d = static_cast<double>(i);
//         double exp_E = E_sum / (i_d + 1.0) / n_cycles_per_check_d / n_spins_d;
//         double exp_Esq = Esq_sum / (i_d + 1.0) / n_cycles_per_check_d / n_spins_d;
//         double exp_ϵ = exp_E / n_spins_d;
//         double exp_ϵsq = exp_Esq / n_spins_d / n_spins_d;
//         double exp_absM = absM_sum / (i_d + 1.0) / n_cycles_per_check_d / n_spins_d;
//         double exp_Msq = Msq_sum / (i_d + 1.0) / n_cycles_per_check_d / n_spins_d;
//         double exp_absm = exp_absM / n_spins_d;
//         double exp_msq = exp_Msq / n_spins_d / n_spins_d;
//         double est_C_V = C_V(n_spins, β, exp_ϵ, exp_ϵsq);
//         double est_χ = χ(n_spins, β, exp_absm, exp_msq);
//     }

//     auto [exp_ϵ_anal, exp_ϵsq_anal] = ϵ_anal(β);
//     auto [exp_absm_anal, exp_msq_anal] = m_anal(β);
//     double est_C_V_anal = C_V_anal(β);
//     double est_χ_anal = χ_anal(β);

//     if (filename != "") {
//         std::ofstream outfile;
//         outfile.open(filename);
//         outfile << "T, ⟨ϵ⟩, ⟨|m|⟩, C_V, χ" << std::endl;

//     }

// }

