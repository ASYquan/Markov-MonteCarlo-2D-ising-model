#include <tuple>
#include <cmath>
#include <cassert>
#include "../include/anal.h"


std::tuple<double,double> ϵ_anal(double β) {
    assert (β > 0.0 && "β must be positive");
    double γ = 3 + std::cosh(8 * β);
    double exp_ϵ_anal = - 2 * std::sinh(8 * β) / γ;
    double exp_ϵsq_anal = 4 * std::cosh(8 * β) / γ;
    return {exp_ϵ_anal, exp_ϵsq_anal};
}

std::tuple<double,double> m_anal(double β) {
    assert (β > 0.0 && "β must be positive");
    double γ = 3 + std::cosh(8 * β);
    double exp_absm_anal = (2 + std::exp(8 * β)) / γ / 2;
    double exp_msq_anal = (1 + std::exp(8 * β)) / γ / 2;
    return {exp_absm_anal, exp_msq_anal};
}

double C_V_anal(double β) {
    assert (β > 0.0 && "β must be positive");
    double γ = 3 + std::cosh(8 * β);
    double exp_C_V_anal = 16 * β * β * (std::cosh(8 * β) / γ - std::sinh(8 * β) *std::sinh(8 * β) / γ / γ);
    return exp_C_V_anal;
}

double χ_anal(double β) {
    assert (β > 0.0 && "β must be positive");
    double γ = 3 + std::cosh(8 * β);
    double exp_χ_anal = β * ((2 + 2 * std::exp(8*β)) / γ - (4 + 4 * std::exp(8 * β) + std::exp(16 * β)) / γ / γ);
    return exp_χ_anal;
}
