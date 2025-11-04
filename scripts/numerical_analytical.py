import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 12, 'ytick.labelsize': 12})

# empty lists to be filled with lists of values
T = []
ϵ = []
m = []
C_V = []
χ = []

# empty lists to be filled with values for each number of cycles
T_n = []
ϵ_n = []
m_n = []
C_V_n = []
χ_n = []

# empty lists to be filled with analytical values
T_a = []
ϵ_a = []
m_a = []
C_V_a = []
χ_a = []

with open("data/ExpectationValues2X2_.txt", "r") as file:
    file.readline()
    lines = file.readlines()
    for line in lines:
        if line.startswith("N"):
            if len(T_n) != 0: 
                # add lists of values to nested lists
                T.append(T_n.copy())
                ϵ.append(ϵ_n.copy())
                m.append(m_n.copy())
                C_V.append(C_V_n.copy())
                χ.append(χ_n.copy())

            T_n = []
            ϵ_n = []
            m_n = []
            C_V_n = []
            χ_n = []
        elif line.startswith("A"):
            # add the last lists of numerical values to nested lists
            T.append(T_n.copy())
            ϵ.append(ϵ_n.copy())
            m.append(m_n.copy())
            C_V.append(C_V_n.copy())
            χ.append(χ_n.copy())

            T_n = []
            ϵ_n = []
            m_n = []
            C_V_n = []
            χ_n = []
        else:
            # add values to n-lists
            values = line.strip().split(",")
            T_n.append(float(values[0]))
            ϵ_n.append(float(values[1]))
            m_n.append(float(values[2]))
            C_V_n.append(float(values[3]))
            χ_n.append(float(values[4]))

    # the last lists to are lists of analytical values
    T_a = T_n.copy()
    ϵ_a = ϵ_n.copy()
    m_a = m_n.copy()
    C_V_a = C_V_n.copy()
    χ_a = χ_n.copy()
        
ϵ_colours = ['#ac0437', '#9d0c43', '#e0115f', '#fa82a7']
m_colours = ['#000f52', '#022f8e', '#1c70c8', '#7dbcde']
C_V_colours = ['#902106', '#d54b00', '#ff8103', '#ffb347']
χ_colours = ['#013220', '#014421', '#8a9a5b', '#bcb88a']
markers = ['|', 'o', 'p', '*']
sizes = [30, 25, 25, 25]

n_cycles = [100, 1000, 10000, 100000]

# plot values
fig = plt.figure(figsize = (12, 8))
gs = fig.add_gridspec(2, 2)
axs = gs.subplots()

for i in range(len(T)): # plot numerical values
    axs[0,0].scatter(T[i], ϵ[i], s=sizes[i], marker=markers[i], color=ϵ_colours[i], label=f'{n_cycles[i]}', zorder=1)
    axs[0,1].scatter(T[i], m[i], s=sizes[i], marker=markers[i], color=m_colours[i], label=f'{n_cycles[i]}', zorder=1)
    axs[1,0].scatter(T[i], C_V[i], s=sizes[i], marker=markers[i], color=C_V_colours[i], label=f'{n_cycles[i]}', zorder=1)
    axs[1,1].scatter(T[i], χ[i], s=sizes[i], marker=markers[i], color=χ_colours[i], label=f'{n_cycles[i]}', zorder=1)

# plot analytical values
axs[0,0].scatter(T_a, ϵ_a, s=30, marker='.', color='k', label=f'Analytical', zorder=1)
axs[0,1].scatter(T_a, m_a, s=30, marker='.', color='k', label=f'Analytical', zorder=1)
axs[1,0].scatter(T_a, C_V_a, s=30, marker='.', color='k', label=f'Analytical', zorder=1)
axs[1,1].scatter(T_a, χ_a, s=30, marker='.', color='k', label=f'Analytical', zorder=1)

axs[0,0].legend(fontsize=12)
axs[0,0].set_ylabel(r'$\left<\epsilon\right>$')
axs[0,0].set_title(r'Energy per spin')

axs[0,1].legend(fontsize=12)
axs[0,1].set_ylabel(r'$\left<|m|\right>$')
axs[0,1].set_title(r'Magnetisation per spin')

axs[1,0].legend(fontsize=12)
axs[1,0].set_ylabel(r'$C_V$')
axs[1,0].set_title(r'Heat capacity')

axs[1,1].legend(fontsize=12)
axs[1,1].set_ylabel(r'$\chi$')
axs[1,1].set_title(r'Susceptibility')

fig.supxlabel(r'Temperature $T$ [$\:J/k_\text{B}$]')
plt.tight_layout()
plt.show()