#include <iostream>
#include <armadillo>
#include <random>
#include <tuple>
#include <chrono>
#include <stdio.h>  
#include <omp.h>
#include "../include/energy.h"
#include "../include/heatcap_suscept.h"
#include "../include/anal.h"
#include "../include/test.h"
#include "../include/MCMC.h"

int main() {
    // Running test before executing main program
    std::cout << "Running tests..." << std::endl;
    std::cout << "-----------------------" << std::endl;
    
    std::cout << "Energy function 1: ";
    test_energy1() ? std::cout << "✅" << std::endl
                   : std::cout << "❌" << std::endl;

    std::cout << "Energy function 2: ";
    test_energy2() ? std::cout << "✅" << std::endl
                   : std::cout << "❌" << std::endl;

    std::cout << "Energy function 3: ";
    test_energy3() ? std::cout << "✅" << std::endl
                   : std::cout << "❌" << std::endl;

    std::cout << "Energy function 4: ";
    test_energy4() ? std::cout << "✅" << std::endl
                   : std::cout << "❌" << std::endl;

    std::cout << "Energy function 5: ";
    test_energy5() ? std::cout << "✅" << std::endl
                   : std::cout << "❌" << std::endl;

    std::cout << "Energy function 6: ";
    test_energy6() ? std::cout << "✅" << std::endl
                   : std::cout << "❌" << std::endl;

    std::cout << "Energy function 7: ";
    test_energy7() ? std::cout << "✅" << std::endl
                   : std::cout << "❌" << std::endl;

    std::cout << "Energy function 8: ";
    test_energy8() ? std::cout << "✅" << std::endl
                   : std::cout << "❌" << std::endl;

    std::cout << "\nTests complete!" << std::endl;
    std::cout << "-----------------------" << std::endl;

    // Creating generator
    // int seed {420'69};
    int seed {1};
    arma::arma_rng::set_seed(seed);
    std::mt19937 generator(seed);

    // Initialising lattice size, number of spins and double-casts
    int L            {}; 
    size_t n_spins   {}; 
    double n_spins_d {};
    double i_d       {};

    // Initialising temporary values
    double E_temp    {}; // E(s)
    double Esq_temp  {}; // (E(s))² 
    double M_temp    {}; // M(s)
    double absM_temp {}; // |M(s)|
    double Msq_temp  {}; // (M(s))² 

    // Initialising sums
    double E_sum     {}; // ΣE(s)
    double Esq_sum   {}; // Σ(E(s))² 
    double absM_sum  {}; // Σ|M(s)|
    double Msq_sum   {}; // Σ(M(s))² 

    // Initialising numerical estimates of expectation values
    double exp_E     {}; // ⟨E⟩
    double exp_Esq   {}; // ⟨E²⟩
    double exp_ϵ     {}; // ⟨ϵ⟩
    double exp_ϵsq   {}; // ⟨ϵ²⟩
    double exp_absM  {}; // ⟨|M|⟩
    double exp_Msq   {}; // ⟨M²⟩
    double exp_absm  {}; // ⟨|m|⟩
    double exp_msq   {}; // ⟨m²⟩
    double est_C_V   {}; // C_V
    double est_χ     {}; // χ

    // Initialising other parameters
    double T         {}; // Temperature
    double β         {}; // Inverse temperature
    double w4        {}; // Probability of spin flip when ΔE = 4J
    double w8        {}; // Probability of spin flip when ΔE = 8J

    // 2X2 lattice case, study numerical results vs. analytical
    {
        // Set lattice size and number of spins
        L = 2; 
        n_spins = L * L; 
        n_spins_d = static_cast<double>(n_spins);

        // Initialising 2x2 lattice
        arma::mat s(L, L); 

        T = 1.0;                   // Temperature
        β = 1.0 / T;               // Inverse temperature
        w4 = std::exp(- β * 4.0);  // Probability of spin flip when ΔE = 4J
        w8 = std::exp(- β * 8.0);  // Probability of spin flip when ΔE = 8J

        // Initialise the number of Monte Carlo cycles to perform and number of times to check expectation values
        size_t n_checks = 10;
        size_t n_cycles_per_check = 1000;
        double n_cycles_per_check_d = static_cast<double>(n_cycles_per_check);

        // Serial
        std::cout << "########################################" << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "Running 2x2 lattice in serial  " << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        auto t1 = std::chrono::high_resolution_clock::now();
        {
            // Fill lattice with values
            s = arma::mat(L, L, arma::fill::randu);   
            s = arma::sign(s - 0.5);

            // Set initial energy and magnetisation values, and start sums
            E_temp = energy(s);   
            Esq_temp = E_temp * E_temp;  
            M_temp = arma::accu(s);
            absM_temp = std::abs(M_temp);   
            Msq_temp = M_temp * M_temp;     
            E_sum = E_temp;    
            Esq_sum = Esq_temp;  
            absM_sum = absM_temp;   
            Msq_sum = Msq_temp;  

            // Perform Monte Carlo steps
            for (size_t i = 0; i < n_checks; ++i) {
                for (size_t j = 0; j < n_cycles_per_check; ++j) {
                    monte_carlo_cycle(s, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, w4, w8, generator); 
                }
                // Calculate expectation values so far
                i_d = static_cast<double>(i);
                exp_E = E_sum / (i_d + 1.0) / n_cycles_per_check_d / n_spins_d;
                exp_Esq = Esq_sum / (i_d + 1.0) / n_cycles_per_check_d / n_spins_d;
                exp_ϵ = exp_E / n_spins_d;
                exp_ϵsq = exp_Esq / n_spins_d / n_spins_d;
                exp_absM = absM_sum / (i_d + 1.0) / n_cycles_per_check_d / n_spins_d;
                exp_Msq = Msq_sum / (i_d + 1.0) / n_cycles_per_check_d / n_spins_d;
                exp_absm = exp_absM / n_spins_d;
                exp_msq = exp_Msq / n_spins_d / n_spins_d;
                est_C_V = C_V(n_spins, β, exp_ϵ, exp_ϵsq);
                est_χ = χ(n_spins, β, exp_absm, exp_msq);

                // Print results
                std::cout << "Results after " << (i + 1) * n_cycles_per_check << " Monte Carlo cycles:" << std::endl;

                printf("⟨E⟩   = %.4f J  ⟨E²⟩ = %.4f J²\n", exp_E, exp_Esq);
                printf("⟨ϵ⟩   = %.4f J  ⟨ϵ²⟩ = %.4f J²\n", exp_ϵ, exp_ϵsq);
                printf("⟨|M|⟩ = %.4f     ⟨M²⟩ = %.4f\n", exp_absM, exp_Msq);
                printf("⟨|m|⟩ = %.4f     ⟨m²⟩ = %.4f\n", exp_absm, exp_msq);
                printf("C_V   = %.4f k_B\n", est_C_V);
                printf("χ     = %.4f J⁻¹\n", est_χ);

                std::cout << "----------------------------------------" << std::endl;

                // RIGHT NOW THERE IS SOMETHING WRONG WITH THE ANALYTICAL SUSCEPTIBILITY VALUE (?) - DOUBLE CHECK THE EXPRESSION
            }
            // Calculate expectation values analytically
            auto [exp_ϵ_anal, exp_ϵsq_anal] = ϵ_anal(β);
            double exp_E_anal = n_spins_d * exp_ϵ_anal;
            double exp_Esq_anal = n_spins_d*n_spins_d * exp_ϵsq_anal;
            auto [exp_absm_anal, exp_msq_anal] = m_anal(β);
            double exp_absM_anal = n_spins_d * exp_absm_anal;
            double exp_Msq_anal = n_spins_d * n_spins_d * exp_msq_anal;
            double est_C_V_anal = C_V_anal(β);
            double est_χ_anal = χ_anal(β);

            std::cout << "----------------------------------------" << std::endl;
            std::cout << "Analytical solution: " << std::endl;

            printf("⟨E⟩   = %.4f J  ⟨E²⟩ = %.4f J²\n", exp_E_anal, exp_Esq_anal);
            printf("⟨ϵ⟩   = %.4f J  ⟨ϵ²⟩ = %.4f J²\n", exp_ϵ_anal, exp_ϵsq_anal);
            printf("⟨|M|⟩ = %.4f     ⟨M²⟩ = %.4f\n", exp_absM_anal, exp_Msq_anal);
            printf("⟨|m|⟩ = %.4f     ⟨m²⟩ = %.4f\n", exp_absm_anal, exp_msq_anal);
            printf("C_V   = %.4f k_B\n", est_C_V_anal);
            printf("χ     = %.4f J⁻¹\n", est_χ_anal);

            std::cout << "----------------------------------------" << std::endl;
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        std::cout << "########################################" << std::endl;
        std::cout << "Time elapsed for 2x2 lattice in serial: " << std::chrono::duration<double>(t2 - t1).count() << " seconds" << std::endl;
        // std::cout << "----------------------------------------" << std::endl;

        

        #ifdef _OPENMP
        // Parallel. Don't print calculations because we're only interested in the time it takes
        // std::cout << "########################################" << std::endl;
        // std::cout << "----------------------------------------" << std::endl;
        // std::cout << "Running 2x2 lattice in parallel" << std::endl;
        // std::cout << "----------------------------------------" << std::endl;
        t1 = std::chrono::high_resolution_clock::now();
        {
            // Fill lattice with values
            s = arma::mat(L, L, arma::fill::randu);   
            s = arma::sign(s - 0.5);

            // Set initial energy and magnetisation values, and start sums
            E_temp = energy(s);   
            Esq_temp = E_temp * E_temp;  
            M_temp = arma::accu(s);
            absM_temp = std::abs(M_temp);   
            Msq_temp = M_temp * M_temp;     
            E_sum = E_temp;    
            Esq_sum = Esq_temp;  
            absM_sum = absM_temp;   
            Msq_sum = Msq_temp;  

            // Perform Monte Carlo steps
            #pragma omp parallel for
            for (size_t i = 0; i < n_checks; ++i) {
                for (size_t j = 0; j < n_cycles_per_check; ++j) {
                    monte_carlo_cycle(s, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, w4, w8, generator); 
                }
                // Calculate expectation values so far
                i_d = static_cast<double>(i);
                exp_E = E_sum / (i_d + 1.0) / n_cycles_per_check_d / n_spins_d;
                exp_Esq = Esq_sum / (i_d + 1.0) / n_cycles_per_check_d / n_spins_d;
                exp_ϵ = exp_E / n_spins_d;
                exp_ϵsq = exp_Esq / n_spins_d / n_spins_d;
                exp_absM = absM_sum / (i_d + 1.0) / n_cycles_per_check_d / n_spins_d;
                exp_Msq = Msq_sum / (i_d + 1.0) / n_cycles_per_check_d / n_spins_d;
                exp_absm = exp_absM / n_spins_d;
                exp_msq = exp_Msq / n_spins_d / n_spins_d;
                est_C_V = C_V(n_spins, β, exp_ϵ, exp_ϵsq);
                est_χ = χ(n_spins, β, exp_absm, exp_msq);
            }
            // Calculate expectation values analytically
            // Only calculated to make runtimes comparable with serial case
            auto [exp_ϵ_anal, exp_ϵsq_anal] = ϵ_anal(β);
            auto [exp_absm_anal, exp_msq_anal] = m_anal(β);
            double est_C_V_anal = C_V_anal(β);
            double est_χ_anal = χ_anal(β);
        }
        t2 = std::chrono::high_resolution_clock::now();
        // std::cout << "----------------------------------------" << std::endl;
        std::cout << "Time elapsed for 2x2 lattice in parallel: " << std::chrono::duration<double>(t2 - t1).count() << " seconds" << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        #endif


        std::ofstream ExpectationValues2X2;
        ExpectationValues2X2.open("data/ExpectationValues2X2.txt");
        ExpectationValues2X2 << "T, ⟨ϵ⟩, ⟨|m|⟩, C_V, χ" << std::endl;

        std::vector<size_t> n_cycles_vals = {100, 1'000, 10'000, 100'000};    // Number of Monte Carlo cycles to perform
        arma::vec T_vals = arma::linspace(1.0, 2.4, 100);                             // Temperature values
        double n_cycles_d {};

        for (auto n_cycles : n_cycles_vals) {
            // std::cout << "########################################" << std::endl;
            // std::cout << "----------------------------------------" << std::endl;
            // std::cout << "Running 2x2 lattice with 100 temperature \nValues in range T ∈ [1.0, 2.4]\n" << n_cycles << " Monte Carlo cycles" << std::endl;
            // std::cout << "----------------------------------------" << std::endl;
            t1 = std::chrono::high_resolution_clock::now();

            n_cycles_d = static_cast<double>(n_cycles);
            ExpectationValues2X2 << "Number of Monte Carlo cycles: " << n_cycles << std::endl;

            for (auto T : T_vals) {
                // Generate LxL lattice with unordered configuration
                s = arma::mat(L, L, arma::fill::randu);   
                s = arma::sign(s - 0.5);

                // Set initial energy and magnetisation values, and start sums
                E_temp = energy(s);   
                Esq_temp = E_temp * E_temp;  
                M_temp = arma::accu(s);
                absM_temp = std::abs(M_temp);   
                Msq_temp = M_temp * M_temp;     
                E_sum = E_temp;    
                Esq_sum = Esq_temp;  
                absM_sum = absM_temp;   
                Msq_sum = Msq_temp;  

                // Set inverse temperature and probabilities
                β = 1.0 / T;
                w4 = std::exp(- β * 4.0);
                w8 = std::exp(- β * 8.0);

                // Perform Monte Carlo steps
                for (size_t i = 0; i < n_cycles; ++i) {
                    monte_carlo_cycle(s, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, w4, w8, generator); 
                }

                // Calculate expectation values
                exp_E = E_sum / n_cycles_d / n_spins_d;
                exp_Esq = Esq_sum / n_cycles_d / n_spins_d;
                exp_ϵ = exp_E / n_spins_d;
                exp_ϵsq = exp_Esq / n_spins_d / n_spins_d;
                exp_absM = absM_sum / n_cycles_d / n_spins_d;
                exp_Msq = Msq_sum / n_cycles_d / n_spins_d;
                exp_absm = exp_absM / n_spins_d;
                exp_msq = exp_Msq / n_spins_d / n_spins_d;
                est_C_V = C_V(n_spins, β, exp_ϵ, exp_ϵsq);
                est_χ = χ(n_spins, β, exp_absm, exp_msq);

                // Write to file
                ExpectationValues2X2 << T << ", " << exp_ϵ << ", " << exp_absm << ", " << est_C_V << ", " << est_χ << std::endl;
            }
            auto t2 = std::chrono::high_resolution_clock::now();
            // std::cout << "----------------------------------------" << std::endl;
            std::cout << "Time elapsed for 2x2 lattice with 100 temperatures ∈ [1.0, 2.4]\n" << n_cycles << " Monte Carlo cycles: " << std::chrono::duration<double>(t2 - t1).count() << " seconds" << std::endl << std::endl;
            // std::cout << "----------------------------------------" << std::endl;
        }
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "Writing analytical solutions to file " << std::endl;

        ExpectationValues2X2 << "Analyticals solution:" << std::endl;

        for (auto T : T_vals) {
            // Set inverse temperature 
            β = 1.0 / T;

            // Calculate expectation values analytically
            auto [exp_ϵ_anal, exp_ϵsq_anal] = ϵ_anal(β);
            auto [exp_absm_anal, exp_msq_anal] = m_anal(β);
            double est_C_V_anal = C_V_anal(β);
            double est_χ_anal = χ_anal(β);

            // Write to file
            ExpectationValues2X2 << T << ", " << exp_ϵ_anal << ", " << exp_absm_anal << ", " << est_C_V_anal << ", " << est_χ_anal << std::endl;
        }
        ExpectationValues2X2.close();
        // std::cout << "----------------------------------------" << std::endl;
        std::cout << "Finished writing analytical solutions to file " << std::endl << std::endl;
        // std::cout << "----------------------------------------" << std::endl;
    }
    // std::exit(1);



    // 20X20 lattice case, study burn-time and probability distribution. 
    // Write ⟨ϵ⟩ and ⟨|m|⟩ to one file for each cycle, and every sampled energy ϵ to another file for each step
    // Do the first for both ordered and unordered, second for only unordered
    // Do this both for T = 1.0 and for T = 2.4
    // Initialise the number of Monte Carlo cycles to perform
    size_t n_cycles = 1'000;
    // TODO maybe change?
    {
        std::cout << "########################################" << std::endl;
        // std::cout << "----------------------------------------" << std::endl;
        std::cout << "Running 20x20 lattice" << std::endl;
        std::cout << "----------------------------------------" << std::endl;

        // Set lattice size and number of spins
        L = 20; 
        n_spins = L * L; 
        n_spins_d = static_cast<double>(n_spins);

        // Initialise 20x20 lattice
        arma::mat s(L, L); 

        // For T = 1.0
        {
            T = 1.0;                  // Temperature
            β = 1.0 / T;              // Inverse temperature
            w4 = std::exp(- β * 4.0); // Probability of spin flip when ΔE = 4J
            w8 = std::exp(- β * 8.0); // Probability of spin flip when ΔE = 8J

            // Ordered configuration
            {
                // Serial
                // std::cout << "########################################" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                // std::cout << "Running T = 1.0 with ordered configuration in serial" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                auto t1 = std::chrono::high_resolution_clock::now();
                {
                    // Generate ordered lattice with all positive spins
                    s = arma::mat(L, L, arma::fill::ones);

                    // Set initial energy and magnetisation values, and start sums
                    E_temp = energy(s);   
                    Esq_temp = E_temp * E_temp;  
                    M_temp = arma::accu(s);
                    absM_temp = std::abs(M_temp);   
                    Msq_temp = M_temp * M_temp;     
                    E_sum = E_temp;    
                    Esq_sum = Esq_temp;  
                    absM_sum = absM_temp;   
                    Msq_sum = Msq_temp;  

                    // Perform Monte Carlo steps and write to file
                    std::ofstream Ordered20X20_T1;
                    Ordered20X20_T1.open("data/Ordered20X20_T1.txt");
                    Ordered20X20_T1 << "⟨ϵ⟩, ⟨|m|⟩" << std::endl;

                    for (size_t i = 0; i < n_cycles; ++i) {
                        monte_carlo_cycle(s, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, w4, w8, generator); 
                        // Calculate expectation values and write to file
                        i_d = static_cast<double>(i);
                        exp_ϵ = E_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        exp_absm = absM_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        Ordered20X20_T1 << exp_ϵ << ", " << exp_absm << std::endl;
                    }
                    Ordered20X20_T1.close();
                }
                auto t2 = std::chrono::high_resolution_clock::now();
                // std::cout << "----------------------------------------" << std::endl;
                std::cout << "Time elapsed for T = 1.0 with ordered configuration in serial: " << std::chrono::duration<double>(t2 - t1).count() << " seconds" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                
                #ifdef _OPENMP
                // Parallel
                // std::cout << "########################################" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                // std::cout << "Running T = 1.0 with ordered configuration in parallel" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                t1 = std::chrono::high_resolution_clock::now();
                {
                    // Generate ordered lattice with all positive spins
                    s = arma::mat(L, L, arma::fill::ones);

                    // Set initial energy and magnetisation values, and start sums
                    E_temp = energy(s);   
                    Esq_temp = E_temp * E_temp;  
                    M_temp = arma::accu(s);
                    absM_temp = std::abs(M_temp);   
                    Msq_temp = M_temp * M_temp;     
                    E_sum = E_temp;    
                    Esq_sum = Esq_temp;  
                    absM_sum = absM_temp;   
                    Msq_sum = Msq_temp;   

                    // Perform Monte Carlo steps and write to file
                    std::ofstream Ordered20X20_T1_par;
                    Ordered20X20_T1_par.open("data/Ordered20X20_T1_par.txt");
                    Ordered20X20_T1_par << "⟨ϵ⟩, ⟨|m|⟩" << std::endl;

                    #pragma omp parallel for private(E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, exp_ϵ, exp_absm)
                    for (size_t i = 0; i < n_cycles; ++i) {
                        monte_carlo_cycle(s, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, w4, w8, generator); 
                        // Calculate expectation values and write to file
                        i_d = static_cast<double>(i);
                        exp_ϵ = E_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        exp_absm = absM_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        Ordered20X20_T1_par << exp_ϵ << ", " << exp_absm << std::endl;
                    }
                    Ordered20X20_T1_par.close();
                }
                t2 = std::chrono::high_resolution_clock::now();
                // std::cout << "----------------------------------------" << std::endl;
                std::cout << "Time elapsed for T = 1.0 with ordered configuration in parallel: " << std::chrono::duration<double>(t2 - t1).count() << " seconds" << std::endl;
                std::cout << "----------------------------------------" << std::endl;
                #endif
            }



            // Unordered configuration and probability function
            {
                // Serial
                // std::cout << "########################################" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                // std::cout << "Running T = 1.0 with unordered configuration and probability function in serial" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                auto t1 = std::chrono::high_resolution_clock::now();
                {
                    // Generate unordered lattice
                    s = arma::mat(L, L, arma::fill::randu);  
                    s = arma::sign(s - 0.5);

                    // Set initial energy and magnetisation values, and start sums
                    E_temp = energy(s);   
                    Esq_temp = E_temp * E_temp;  
                    M_temp = arma::accu(s);
                    absM_temp = std::abs(M_temp);   
                    Msq_temp = M_temp * M_temp;     
                    E_sum = E_temp;    
                    Esq_sum = Esq_temp;  
                    absM_sum = absM_temp;   
                    Msq_sum = Msq_temp;   

                    // Create file containing ϵ-values to be used for estimating probability distribution
                    std::ofstream Probability20X20ϵ_T1;
                    Probability20X20ϵ_T1.open("data/Probability20X20ϵ_T1.txt");
                    Probability20X20ϵ_T1 << "ϵ" << std::endl;
                    Probability20X20ϵ_T1.close();

                    // Create file containing m-values to be used for studying discrepancies
                    std::ofstream Probability20X20m_T1;
                    Probability20X20m_T1.open("data/Probability20X20m_T1.txt");
                    Probability20X20m_T1 << "m" << std::endl;
                    Probability20X20m_T1.close();

                    // Perform Monte Carlo steps and write to file
                    std::ofstream Unordered20X20_T1;
                    Unordered20X20_T1.open("data/Unordered20X20_T1.txt");
                    Unordered20X20_T1 << "⟨ϵ⟩, ⟨|m|⟩" << std::endl;

                    for (size_t i = 0; i < n_cycles; ++i) {
                        monte_carlo_cycle(s, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, w4, w8, generator, "data/Probability20X20ϵ_T1.txt", "data/Probability20X20m_T1.txt"); 
                        // Calculate expectation values and write to file
                        i_d = static_cast<double>(i);
                        exp_ϵ = E_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        exp_absm = absM_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        Unordered20X20_T1 << exp_ϵ << ", " << exp_absm << std::endl;
                    }
                    Unordered20X20_T1.close();
                }
                auto t2 = std::chrono::high_resolution_clock::now();
                // std::cout << "----------------------------------------" << std::endl;
                std::cout << "Time elapsed for T = 1.0 with unordered configuration and probability function in serial: " << std::chrono::duration<double>(t2 - t1).count() << " seconds" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;

                #ifdef _OPENMP
                // Parallel
                // std::cout << "########################################" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                // std::cout << "Running T = 1.0 with unordered configuration and probability function in parallel" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                t1 = std::chrono::high_resolution_clock::now();
                {
                    // Generate unordered lattice
                    s = arma::mat(L, L, arma::fill::randu);  
                    s = arma::sign(s - 0.5);

                    // Set initial energy and magnetisation values, and start sums
                    E_temp = energy(s);   
                    Esq_temp = E_temp * E_temp;  
                    M_temp = arma::accu(s);
                    absM_temp = std::abs(M_temp);   
                    Msq_temp = M_temp * M_temp;     
                    E_sum = E_temp;    
                    Esq_sum = Esq_temp;  
                    absM_sum = absM_temp;   
                    Msq_sum = Msq_temp;  

                    // Create file containing ϵ-values to be used for estimating probability distribution
                    std::ofstream Probability20X20ϵ_T1_par;
                    Probability20X20ϵ_T1_par.open("data/Probability20X20ϵ_T1_par.txt");
                    Probability20X20ϵ_T1_par << "ϵ" << std::endl;
                    Probability20X20ϵ_T1_par.close();

                    // Create file containing m-values to be used for studying discrepancies
                    std::ofstream Probability20X20m_T1_par;
                    Probability20X20m_T1_par.open("data/Probability20X20m_T1_par.txt");
                    Probability20X20m_T1_par << "m" << std::endl;
                    Probability20X20m_T1_par.close();

                    // Perform Monte Carlo steps and write to file
                    std::ofstream Unordered20X20_T1_par;
                    Unordered20X20_T1_par.open("data/Unordered20X20_T1_par.txt");
                    Unordered20X20_T1_par << "⟨ϵ⟩, ⟨|m|⟩" << std::endl;

                    #pragma omp parallel for private(E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, exp_ϵ, exp_absm)
                    for (size_t i = 0; i < n_cycles; ++i) {
                        monte_carlo_cycle(s, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, w4, w8, generator, "data/Probability20X20ϵ_T1_par.txt", "data/Probability20X20m_T1_par.txt"); 
                        // Calculate expectation values and write to file
                        i_d = static_cast<double>(i);
                        exp_ϵ = E_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        exp_absm = absM_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        Unordered20X20_T1_par << exp_ϵ << ", " << exp_absm << std::endl;
                    }
                    Unordered20X20_T1_par.close();
                }
                t2 = std::chrono::high_resolution_clock::now();
                // std::cout << "----------------------------------------" << std::endl;
                std::cout << "Time elapsed for T = 1.0 with unordered configuration and probability function in parallel: " << std::chrono::duration<double>(t2 - t1).count() << " seconds" << std::endl;
                std::cout << "----------------------------------------" << std::endl;
                #endif
            }
        }





        // For T = 2.4
        {
            T = 2.4;                  // Temperature
            β = 1.0 / T;              // Inverse temperature
            w4 = std::exp(- β * 4.0); // Probability of spin flip when ΔE = 4J
            w8 = std::exp(- β * 8.0); // Probability of spin flip when ΔE = 8J

            // Ordered configuration
            {
                // Serial
                // std::cout << "########################################" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                // std::cout << "Running T = 2.4 with ordered configuration in serial" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                auto t1 = std::chrono::high_resolution_clock::now();
                {
                    // Generate ordered lattice with all positive spins
                    s = arma::mat(L, L, arma::fill::ones); 

                    // Set initial energy and magnetisation values, and start sums
                    E_temp = energy(s);   
                    Esq_temp = E_temp * E_temp;  
                    M_temp = arma::accu(s);
                    absM_temp = std::abs(M_temp);   
                    Msq_temp = M_temp * M_temp;     
                    E_sum = E_temp;    
                    Esq_sum = Esq_temp;  
                    absM_sum = absM_temp;   
                    Msq_sum = Msq_temp;   

                    // Perform Monte Carlo steps and write to file
                    std::ofstream Ordered20X20_T2;
                    Ordered20X20_T2.open("data/Ordered20X20_T2.txt");
                    Ordered20X20_T2 << "⟨ϵ⟩, ⟨|m|⟩" << std::endl;

                    for (size_t i = 0; i < n_cycles; ++i) {
                        monte_carlo_cycle(s, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, w4, w8, generator); 
                        // Calculate expectation values and write to file
                        i_d = static_cast<double>(i);
                        exp_ϵ = E_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        exp_absm = absM_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        Ordered20X20_T2 << exp_ϵ << ", " << exp_absm << std::endl;
                    }
                    Ordered20X20_T2.close();
                }
                auto t2 = std::chrono::high_resolution_clock::now();
                // std::cout << "----------------------------------------" << std::endl;
                std::cout << "Time elapsed for T = 2.4 with ordered configuration in serial: " << std::chrono::duration<double>(t2 - t1).count() << " seconds" << std::endl;

                #ifdef _OPENMP
                // Parallel
                // std::cout << "########################################" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                // std::cout << "Running T = 2.4 with ordered configuration in parallel" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                t1 = std::chrono::high_resolution_clock::now();
                {
                    // Generate ordered lattice with all positive spins
                    s = arma::mat(L, L, arma::fill::ones);  

                    // Set initial energy and magnetisation values, and start sums
                    E_temp = energy(s);   
                    Esq_temp = E_temp * E_temp;  
                    M_temp = arma::accu(s);
                    absM_temp = std::abs(M_temp);   
                    Msq_temp = M_temp * M_temp;     
                    E_sum = E_temp;    
                    Esq_sum = Esq_temp;  
                    absM_sum = absM_temp;   
                    Msq_sum = Msq_temp;  

                    // Perform Monte Carlo steps and write to file
                    std::ofstream Ordered20X20_T2_par;
                    Ordered20X20_T2_par.open("data/Ordered20X20_T2_par.txt");
                    Ordered20X20_T2_par << "⟨ϵ⟩, ⟨|m|⟩" << std::endl;

                    #pragma omp parallel for private(E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, exp_ϵ, exp_absm)
                    for (size_t i = 0; i < n_cycles; ++i) {
                        monte_carlo_cycle(s, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, w4, w8, generator); 
                        // Calculate expectation values and write to file
                        i_d = static_cast<double>(i);
                        exp_ϵ = E_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        exp_absm = absM_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        Ordered20X20_T2_par << exp_ϵ << ", " << exp_absm << std::endl;
                    }
                    Ordered20X20_T2_par.close();
                }
                t2 = std::chrono::high_resolution_clock::now();
                // std::cout << "----------------------------------------" << std::endl;
                std::cout << "Time elapsed for T = 2.4 with ordered configuration in parallel: " << std::chrono::duration<double>(t2 - t1).count() << " seconds" << std::endl;
                std::cout << "----------------------------------------" << std::endl;
                #endif
            }



            // Unordered configuration and probability function
            {
                // Serial
                // std::cout << "########################################" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                // std::cout << "Running T = 2.4 with unordered configuration and probability function in serial" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                auto t1 = std::chrono::high_resolution_clock::now();
                {
                    // Generate unordered lattice
                    s = arma::mat(L, L, arma::fill::randu);  
                    s = arma::sign(s - 0.5);

                    // Set initial energy and magnetisation values, and start sums
                    E_temp = energy(s);   
                    Esq_temp = E_temp * E_temp;  
                    M_temp = arma::accu(s);
                    absM_temp = std::abs(M_temp);   
                    Msq_temp = M_temp * M_temp;     
                    E_sum = E_temp;    
                    Esq_sum = Esq_temp;  
                    absM_sum = absM_temp;   
                    Msq_sum = Msq_temp;  

                    // Create file containing ϵ-values to be used for estimating probability distribution
                    std::ofstream Probability20X20ϵ_T2;
                    Probability20X20ϵ_T2.open("data/Probability20X20ϵ_T2.txt");
                    Probability20X20ϵ_T2 << "ϵ" << std::endl;
                    Probability20X20ϵ_T2.close();

                    // Create file containing m-values to be used for studying discrepancies
                    std::ofstream Probability20X20m_T2;
                    Probability20X20m_T2.open("data/Probability20X20m_T2.txt");
                    Probability20X20m_T2 << "m" << std::endl;
                    Probability20X20m_T2.close();

                    // Perform Monte Carlo steps and write to file
                    std::ofstream Unordered20X20_T2;
                    Unordered20X20_T2.open("data/Unordered20X20_T2.txt");
                    Unordered20X20_T2 << "⟨ϵ⟩, ⟨|m|⟩" << std::endl;

                    for (size_t i = 0; i < n_cycles; ++i) {
                        monte_carlo_cycle(s, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, w4, w8, generator, "data/Probability20X20ϵ_T2.txt", "data/Probability20X20m_T2.txt"); 
                        // Calculate expectation values and write to file
                        i_d = static_cast<double>(i);
                        exp_ϵ = E_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        exp_absm = absM_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        Unordered20X20_T2 << exp_ϵ << ", " << exp_absm << std::endl;
                    }
                    Unordered20X20_T2.close();
                }
                auto t2 = std::chrono::high_resolution_clock::now();
                // std::cout << "----------------------------------------" << std::endl;
                std::cout << "Time elapsed for T = 2.4 with unordered configuration and probability function in serial: " << std::chrono::duration<double>(t2 - t1).count() << " seconds" << std::endl;
                // std::cout << "----------------------------------------" << std::endl << std::endl;

                #ifdef _OPENMP
                // Parallel
                // std::cout << "########################################" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                // std::cout << "Running T = 2.4 with unordered configuration and probability function in parallel" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                t1 = std::chrono::high_resolution_clock::now();
                {
                    // Generate unordered lattice
                    s = arma::mat(L, L, arma::fill::randu);  
                    s = arma::sign(s - 0.5);

                    // Set initial energy and magnetisation values, and start sums
                    E_temp = energy(s);   
                    Esq_temp = E_temp * E_temp;  
                    M_temp = arma::accu(s);
                    absM_temp = std::abs(M_temp);   
                    Msq_temp = M_temp * M_temp;     
                    E_sum = E_temp;    
                    Esq_sum = Esq_temp;  
                    absM_sum = absM_temp;   
                    Msq_sum = Msq_temp;   

                    // Create file containing ϵ-values to be used for estimating probability distribution
                    std::ofstream Probability20X20ϵ_T2_par;
                    Probability20X20ϵ_T2_par.open("data/Probability20X20ϵ_T2_par.txt");
                    Probability20X20ϵ_T2_par << "ϵ" << std::endl;
                    Probability20X20ϵ_T2_par.close();

                    // Create file containing m-values to be used for studying discrepancies
                    std::ofstream Probability20X20m_T2_par;
                    Probability20X20m_T2_par.open("data/Probability20X20m_T2_par.txt");
                    Probability20X20m_T2_par << "m" << std::endl;
                    Probability20X20m_T2_par.close();

                    // Perform Monte Carlo steps and write to file
                    std::ofstream Unordered20X20_T2_par;
                    Unordered20X20_T2_par.open("data/Unordered20X20_T2_par.txt");
                    Unordered20X20_T2_par << "⟨ϵ⟩, ⟨|m|⟩" << std::endl;
                    
                    #pragma omp parallel for private(E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, exp_ϵ, exp_absm)
                    for (size_t i = 0; i < n_cycles; ++i) {
                        monte_carlo_cycle(s, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, w4, w8, generator, "data/Probability20X20ϵ_T2_par.txt", "data/Probability20X20m_T2_par.txt"); 
                        // Calculate expectation values and write to file
                        i_d = static_cast<double>(i);
                        exp_ϵ = E_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        exp_absm = absM_sum / (i_d + 1.0) / n_spins_d / n_spins_d;
                        Unordered20X20_T2_par << exp_ϵ << ", " << exp_absm << std::endl;
                    }
                    Unordered20X20_T2_par.close();
                }
                t2 = std::chrono::high_resolution_clock::now();
                // std::cout << "----------------------------------------" << std::endl;
                std::cout << "Time elapsed for T = 2.4 with unordered configuration and probability function in parallel: " << std::chrono::duration<double>(t2 - t1).count() << " seconds" << std::endl;
                std::cout << "----------------------------------------" << std::endl << std::endl;
                #endif
            }
        }
    }



    std::cout << "########################################" << std::endl;
    std::cout << "########################################" << std::endl;
    std::cout << "########################################" << std::endl;

    // Phase transitions
    {
        std::vector<size_t> L_vals = {40, 60, 80, 100};        // Lattice sizes
        arma::vec T_vals = arma::linspace(2.1, 2.4, 20);    // Temperature values
        
        size_t n_cycles = 50'000;    // Number of Monte Carlo cycles to perform

        #ifdef _OPENMP
        {

            for (auto L : L_vals) {
                // std::cout << "########################################" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                // std::cout << "Running " << L << "x" << L << " lattice with 20 temperature values in range T ∈ [2.1, 2.4]" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                
                auto t1 = std::chrono::high_resolution_clock::now();


                // Set number of spins
                n_spins = L * L; 
                n_spins_d = static_cast<double>(n_spins);

                arma::mat s(L, L);
                #pragma omp parallel for private(s, β, w4, w8, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, exp_E, exp_Esq, exp_ϵ, exp_absM, exp_Msq, exp_absm, est_C_V, est_χ) num_threads(20)
                for (auto T : T_vals) {
                    // Generate LxL lattice with unordered configuration
                    s = arma::mat(L, L, arma::fill::randu);   
                    s = arma::sign(s - 0.5);

                    // Set initial energy and magnetisation values, and start sums
                    E_temp = energy(s);   
                    Esq_temp = E_temp * E_temp;  
                    M_temp = arma::accu(s);
                    absM_temp = std::abs(M_temp);   
                    Msq_temp = M_temp * M_temp;     
                    E_sum = E_temp;    
                    Esq_sum = Esq_temp;  
                    absM_sum = absM_temp;   
                    Msq_sum = Msq_temp;  

                    // Set inverse temperature and probabilities
                    β = 1.0 / T;
                    w4 = std::exp(- β * 4.0);
                    w8 = std::exp(- β * 8.0);

                    // Perform Monte Carlo steps
                    for (size_t i = 0; i < n_cycles; ++i) {
                        monte_carlo_cycle(s, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, w4, w8, generator); 
                    }

                    // Calculate expectation values
                    double n_cycles_d =  static_cast<double>(n_cycles);
                    exp_E = E_sum / n_cycles_d / n_spins_d;
                    exp_Esq = Esq_sum / n_cycles_d / n_spins_d;
                    exp_ϵ = exp_E / n_spins_d;
                    exp_ϵsq = exp_Esq / n_spins_d / n_spins_d;
                    exp_absM = absM_sum / n_cycles_d / n_spins_d;
                    exp_Msq = Msq_sum / n_cycles_d / n_spins_d;
                    exp_absm = exp_absM / n_spins_d;
                    exp_msq = exp_Msq / n_spins_d / n_spins_d;
                    est_C_V = C_V(n_spins, β, exp_ϵ, exp_ϵsq);
                    est_χ = χ(n_spins, β, exp_absm, exp_msq);

                    // Write to file


                    //TODO write to file outside loop instead?
                }
                auto t2 = std::chrono::high_resolution_clock::now();
                // std::cout << "----------------------------------------" << std::endl;
                std::cout << "Parallell: Time elapsed for " << L << "x" << L << " lattice with 20 temperature values ∈ [2.1, 2.4]: " << std::chrono::duration<double>(t2 - t1).count() << " seconds" << std::endl;
                }   
                std::cout << "----------------------------------------" << std::endl;
        }
        #endif



        {
            T_vals = arma::linspace(2.20, 2.35, 20); 

            for (auto L : L_vals) {
                // std::cout << "########################################" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                // std::cout << "Running " << L << "x" << L << " lattice with 20 temperature values in range T ∈ [2.1, 2.4]" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                
                auto t1 = std::chrono::high_resolution_clock::now();

                // Set number of spins
                n_spins = L * L; 
                n_spins_d = static_cast<double>(n_spins);

                arma::mat s(L, L);
                for (auto T : T_vals) {
                    // Generate LxL lattice with unordered configuration
                    s = arma::mat(L, L, arma::fill::randu);   
                    s = arma::sign(s - 0.5);

                    // Set initial energy and magnetisation values, and start sums
                    E_temp = energy(s);   
                    Esq_temp = E_temp * E_temp;  
                    M_temp = arma::accu(s);
                    absM_temp = std::abs(M_temp);   
                    Msq_temp = M_temp * M_temp;     
                    E_sum = E_temp;    
                    Esq_sum = Esq_temp;  
                    absM_sum = absM_temp;   
                    Msq_sum = Msq_temp;  

                    // Set inverse temperature and probabilities
                    β = 1.0 / T;
                    w4 = std::exp(- β * 4.0);
                    w8 = std::exp(- β * 8.0);

                    // Perform Monte Carlo steps
                    for (size_t i = 0; i < n_cycles; ++i) {
                        monte_carlo_cycle(s, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, w4, w8, generator); 
                    }

                    // Calculate expectation values
                    double n_cycles_d =  static_cast<double>(n_cycles);
                    exp_E = E_sum / n_cycles_d / n_spins_d;
                    exp_Esq = Esq_sum / n_cycles_d / n_spins_d;
                    exp_ϵ = exp_E / n_spins_d;
                    exp_ϵsq = exp_Esq / n_spins_d / n_spins_d;
                    exp_absM = absM_sum / n_cycles_d / n_spins_d;
                    exp_Msq = Msq_sum / n_cycles_d / n_spins_d;
                    exp_absm = exp_absM / n_spins_d;
                    exp_msq = exp_Msq / n_spins_d / n_spins_d;
                    est_C_V = C_V(n_spins, β, exp_ϵ, exp_ϵsq);
                    est_χ = χ(n_spins, β, exp_absm, exp_msq);

                    // Write to file

                    //TODO write to file outside loop instead?
                }
                auto t2 = std::chrono::high_resolution_clock::now();
                // std::cout << "----------------------------------------" << std::endl;
                std::cout << "Serial: Time elapsed for " << L << "x" << L << " lattice with 20 temperature values ∈ [2.1, 2.4]: " << std::chrono::duration<double>(t2 - t1).count() << " seconds" << std::endl;
                }   
                std::cout << "----------------------------------------" << std::endl;
        }











        
        // Studying a narrower temperature range
        T_vals = arma::linspace(2.20, 2.35, 50);    // Temperature values

        #ifdef _OPENMP
        {
            std::ofstream ExpectationValuesLXL_Narrow;
            ExpectationValuesLXL_Narrow.open("data/ExpectationValuesLXL_Narrow.txt");
            ExpectationValuesLXL_Narrow << "T, ⟨ϵ⟩, ⟨|m|⟩, C_V, χ" << std::endl;

            for (auto L : L_vals) {
                // std::cout << "########################################" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                // std::cout << "Running " << L << "x" << L << " lattice with 20 temperature values in range T ∈ [2.25, 2.35]" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
                auto t1 = std::chrono::high_resolution_clock::now();

                ExpectationValuesLXL_Narrow << "Lattice size: " << L << "×" << L << std::endl;

                // Set number of spins
                n_spins = L * L; 
                n_spins_d = static_cast<double>(n_spins);

                arma::mat s(L, L);
                #pragma omp parallel for private(s, β, w4, w8, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, exp_E, exp_Esq, exp_ϵ, exp_absM, exp_Msq, exp_absm, est_C_V, est_χ) num_threads(20)
                for (auto T : T_vals) {
                    // Generate LXL lattice with unordered configuration
                    s = arma::mat(L, L, arma::fill::randu);   
                    s = arma::sign(s - 0.5);

                    // Set initial energy and magnetisation values, and start sums
                    E_temp = energy(s);   
                    Esq_temp = E_temp * E_temp;  
                    M_temp = arma::accu(s);
                    absM_temp = std::abs(M_temp);   
                    Msq_temp = M_temp * M_temp;     
                    E_sum = E_temp;    
                    Esq_sum = Esq_temp;  
                    absM_sum = absM_temp;   
                    Msq_sum = Msq_temp;   

                    // Set inverse temperature and probabilities
                    β = 1.0 / T;
                    w4 = std::exp(- β * 4.0);
                    w8 = std::exp(- β * 8.0);

                    // Perform Monte Carlo steps
                    for (size_t i = 0; i < n_cycles; ++i) {
                        monte_carlo_cycle(s, E_temp, Esq_temp, M_temp, absM_temp, Msq_temp, E_sum, Esq_sum, absM_sum, Msq_sum, w4, w8, generator);
                    }

                    // Calculate expectation values
                    double n_cycles_d = static_cast<double>(n_cycles);
                    exp_E = E_sum / n_cycles_d / n_spins_d;
                    exp_Esq = Esq_sum / n_cycles_d / n_spins_d;
                    exp_ϵ = exp_E / n_spins_d;
                    exp_ϵsq = exp_Esq / n_spins_d / n_spins_d;
                    exp_absM = absM_sum / n_cycles_d / n_spins_d;
                    exp_Msq = Msq_sum / n_cycles_d / n_spins_d;
                    exp_absm = exp_absM / n_spins_d;
                    exp_msq = exp_Msq / n_spins_d / n_spins_d;
                    est_C_V = C_V(n_spins, β, exp_ϵ, exp_ϵsq);
                    est_χ = χ(n_spins, β, exp_absm, exp_msq);

                    // Write to file
                    ExpectationValuesLXL_Narrow << T << ", " << exp_ϵ << ", " << exp_absm << ", " << est_C_V << ", " << est_χ << std::endl;
                    //TODO write to file outside loop instead?
                }
                auto t2 = std::chrono::high_resolution_clock::now();
                // std::cout << "----------------------------------------" << std::endl;
                std::cout << "Time elapsed for " << L << "x" << L << " lattice with 50 temperature values ∈ [2.25, 2.35]: " << std::chrono::duration<double>(t2 - t1).count() << " seconds" << std::endl;
                // std::cout << "----------------------------------------" << std::endl;
            }
            ExpectationValuesLXL_Narrow.close();
        }
        #endif
    }
    
    return 0;
}