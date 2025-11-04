# 2D Ising Model — Simple Monte Carlo, Clear Phase Transitions

Curiosity project exploring a classic lattice model of magnetism with a simple sampler and clean plots. 

NB: The data was too large.... had to seperate into folders. 
This repo is for archival purposes for code in Python, LaTex and C++

## What this does (short)
- Simulates a 2D grid of spins (each spin is either +1 or −1).
- Uses Markov Chain Monte Carlo (Metropolis) to sample configurations at a given temperature.
- Tracks energy, magnetization, heat capacity, and susceptibility.
- Detects the phase transition by scanning temperature and estimating the critical temperature.
- Optional parallelization to speed up long runs.

This model can be used in finance, because it is a neat proxy for “regime changes” — at low temperature spins align (ordered), at high temperature they randomize (disordered). That mirrors how markets can flip between calm and turbulent regimes. MCMC sampling is also a staple in Bayesian portfolio models and risk estimation.
Furthermore, MCMC is a workhorse for sampling complex, high‑dimensional distributions. The Ising model gives a clean sandbox to practice sampling, autocorrelation, burn‑in, and parallel runs. Skills that transfer to probabilistic modeling and generative methods.

## What’s included (typical structure)
- Lattice setup with periodic boundaries (every spin has 4 neighbors).
- Metropolis updates: propose a flip, accept with the Boltzmann factor.
- Observables:
  - Energy per spin and magnetization per spin.
  - Heat capacity and susceptibility (from fluctuations).
- Experiments:
  - Burn‑in (equilibration) plots from ordered vs random starts.
  - Histograms of magnetization for different temperatures.
  - Temperature scans across the phase transition.
  - Finite‑size scaling to estimate the critical temperature.


## Typical workflow
- Choose lattice size (e.g., L = 20, 40, 60) and temperature range.
- Run MCMC:
  - Burn‑in until observables stabilize.
  - Collect samples to compute averages and fluctuations.
- Plot:
  - Energy and magnetization vs Monte Carlo cycles (burn‑in).
  - Magnetization histograms at selected temperatures.
  - Energy, |magnetization|, heat capacity, susceptibility vs temperature.
- Estimate the critical temperature via finite‑size scaling.

## Notes
- Use precomputed Boltzmann factors for speed (only five possible energy changes for single‑spin flips).
- Periodic boundaries can be implemented efficiently without if‑statements (index wrapping).
- Parallel speed‑ups: run multiple temperatures or multiple independent chains simultaneously.


## Build and Run
### Running the Executable
From the root directory:
```bash
./src/main
```
### Building with make
To build the project, run `make` in the root folder

```bash
make
```

### Building Without make in the Terminal
#### As a Single Command
```bash
g++ --std=c++17 -Wall -O3 -c src/*.cpp && mv *.o objectFiles && g++ --std=c++17 -Wall -O3 -o src/main objectFiles/*.o && ./src/main
```


#### Each individual Command
1. Compile
```bash
g++ --std=c++17 -Wall -O3 -c src/*.cpp
```
2. Move object files to objectFiles folder
```bash
mv *.o objectFiles
```
3. Link object files and create executable
```bash
g++ --std=c++17 -Wall -O3 -o src/main objectFiles/*.o
```	
4. Run executable
```bash
./src/main
```


