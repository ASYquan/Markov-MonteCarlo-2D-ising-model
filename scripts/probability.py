import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 12, 'ytick.labelsize': 12})

T = [1.0, 2.4]

ϵ_colours = ['tomato', 'mediumvioletred']
m_colours = ['mediumpurple', 'seagreen']

fig = plt.figure(figsize = (12, 8))
gs = fig.add_gridspec(2, 2)
axs = gs.subplots()

for i in range(2):
    ϵ = np.loadtxt(f"data/Probability20X20_T{i+1}.txt", skiprows=1)
    m = np.loadtxt(f"data/Probability20X20m_T{i+1}.txt", skiprows=1)

    Boltzmann = np.exp(-ϵ/T[i])       # Boltzmann factors
    Z = np.sum(Boltzmann)             # partition function

    # create normalized histogram for the current temperature p = Boltzmann/Z
    p = Boltzmann/Z

    steps = np.arange(len(p))
    num_bins_ϵ = len(set(ϵ))
    num_bins_m = len(set(m))

    # create a histogram-like bar plot for current temperature
    axs[0, i].hist(ϵ, bins=num_bins_ϵ, weights=p, color=ϵ_colours[i])
    axs[0, i].set_title(r'$T=$' + f'{T[i]}' + r'$\:J/k_\text{B}$')
    axs[1, i].hist(m, bins=num_bins_m, weights=p, color=m_colours[i])

axs[0, 0].set_xlabel(95*' ' + r'Energy per spin $\varepsilon$ [$\:J$]')
axs[1, 0].set_xlabel(95*' ' + r'Magnetisation per spin $m$')
axs[0, 0].set_ylabel(r'$p_\varepsilon(\varepsilon;T)$')
axs[1, 0].set_ylabel(r'$p_{m}(m;T)$')

plt.tight_layout()
plt.show()