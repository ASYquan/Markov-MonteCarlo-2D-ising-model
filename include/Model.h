#include <cstddef>
#include <string>

class Model {
private:
    size_t L;
    double T;
    size_t n_checks;
    size_t n_cycles_per_check;
    int seed;
public:
    Model(std::size_t L, double T, std::size_t n_checks, std::size_t n_cycles_per_check, int seed);
    ~Model();

    void MCMC(bool parallel = false, bool ordered = false, std::string filename = "");
};