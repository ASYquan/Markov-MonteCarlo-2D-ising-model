#include <cassert>
#include "../include/heatcap_suscept.h"

double C_V(std::size_t n_spins, double β, double exp_ϵ, double exp_ϵsq) {
    assert (β > 0.0 && "β must be positive");
    double n_spins_d = static_cast<double>(n_spins);
    return n_spins_d * β * β * (exp_ϵsq - exp_ϵ * exp_ϵ);
}

double χ(std::size_t n_spins, double β, double exp_absm, double exp_msq) {
    assert (β > 0.0 && "β must be positive");
    double n_spins_d = static_cast<double>(n_spins);
    return n_spins_d * β * (exp_msq - exp_absm * exp_absm);
}