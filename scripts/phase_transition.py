import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 12, 'ytick.labelsize': 12})

# empty lists to be filled with lists of values
T = []
ϵ = []
m = []
C_V = []
χ = []

# empty lists to be filled with values for each lattice
T_L = []
ϵ_L = []
m_L = []
C_V_L = []
χ_L = []

with open("data/ExpectationValuesLXL_Narrow.txt", "r") as file:
#with open("data/ExpectationValuesLXL.txt", "r") as file:
    file.readline()
    lines = file.readlines()
    for line in lines:
        if line.startswith("L"):
            if len(T_L) != 0: 
                # add lists of values to nested lists
                T.append(T_L.copy())
                ϵ.append(ϵ_L.copy())
                m.append(m_L.copy())
                C_V.append(C_V_L.copy())
                χ.append(χ_L.copy())
            T_L = []
            ϵ_L = []
            m_L = []
            C_V_L = []
            χ_L = []
        else:
            # add values to L-lists
            values = line.strip().split(",")
            T_L.append(float(values[0]))
            ϵ_L.append(float(values[1]))
            m_L.append(float(values[2]))
            C_V_L.append(float(values[3]))
            χ_L.append(float(values[4]))
    # add the last lists to nested lists
    T.append(T_L.copy())
    ϵ.append(ϵ_L.copy())
    m.append(m_L.copy())
    C_V.append(C_V_L.copy())
    χ.append(χ_L.copy())

T_sorted = np.sort(T) # lists of temperatures sorted from smallest to largest
ϵ_sorted = []
m_sorted = []
C_V_sorted = []
χ_sorted = []

# sort nested lists of values
for i in range(len(T_sorted)):
    for j in range(len(T_sorted[i])):
        idx = T[i].index(T_sorted[i,j])
        ϵ_L[j] = ϵ[i][idx]
        m_L[j] = m[i][idx]
        C_V_L[j] = C_V[i][idx]
        χ_L[j] = χ[i][idx]

    ϵ_sorted.append(ϵ_L.copy())
    m_sorted.append(m_L.copy())
    C_V_sorted.append(C_V_L.copy())
    χ_sorted.append(χ_L.copy())

# replace old nested lists with the sorted ones
T = np.array(T_sorted)
ϵ = np.array(ϵ_sorted)
m = np.array(m_sorted)
C_V = np.array(C_V_sorted)
χ = np.array(χ_sorted)
        
ϵ_colours = ['#ac0437', '#9d0c43', '#e0115f', '#fa82a7']
m_colours = ['#000f52', '#022f8e', '#1c70c8', '#7dbcde']
C_V_colours = ['#902106', '#d54b00', '#ff8103', '#ffb347']
χ_colours = ['#013220', '#014421', '#8a9a5b', '#bcb88a']
markers = ['x', 'o', 'p', '*']

L = [40, 60, 80, 100]

# plot values
fig = plt.figure(figsize = (12, 8))
gs = fig.add_gridspec(2, 2)
axs = gs.subplots()

for i in range(len(T)):
    axs[0,0].scatter(T[i], ϵ[i], s=30, marker=markers[i], color=ϵ_colours[i], label=f'{L[i]}×{L[i]}', zorder=1)
    axs[0,0].plot(T[i], ϵ[i], color=ϵ_colours[i], alpha=0.3, zorder=0)

    axs[0,1].scatter(T[i], m[i], s=30, marker=markers[i], color=m_colours[i], label=f'{L[i]}×{L[i]}', zorder=1)
    axs[0,1].plot(T[i], m[i], color=m_colours[i], alpha=0.3, zorder=0)

    axs[1,0].scatter(T[i], C_V[i], s=30, marker=markers[i], color=C_V_colours[i], label=f'{L[i]}×{L[i]}', zorder=1)
    axs[1,0].plot(T[i], C_V[i], color=C_V_colours[i], alpha=0.3, zorder=0)

    axs[1,1].scatter(T[i], χ[i], s=30, marker=markers[i], color=χ_colours[i], label=f'{L[i]}×{L[i]}', zorder=1)
    axs[1,1].plot(T[i], χ[i], color=χ_colours[i], alpha=0.3, zorder=0)

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




''' Estimating the critical temperatures of an infinite 2D Ising model '''

T_c = []

for i in range(len(T)):
    # estimate critical temperature of LxL lattice from C_V and χ values
    idx1 = np.argmax(C_V[i])
    idx2 = np.argmax(χ[i])
    T_c.append((T[i, idx1] + T[i, idx2]) / 2) 

L = np.array(L)
T_c = np.array(T_c)
print(T_c)

T_c_inf_est = T_c - 1 / L

T_c_inf = np.mean(T_c_inf_est)      # rough estimate of critical temperature for infinite lattice
T_c_inf_err = np.std(T_c_inf_est)   # rough estimate of uncertainty in critical temperature for infinite lattice
print(T_c_inf, T_c_inf_err)
# 2.258958333333333 0.003420721780560296