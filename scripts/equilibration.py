import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size': 14, 'xtick.labelsize': 12, 'ytick.labelsize': 12})

# filenames of files to extract data from
filenames = ["data/Ordered20X20_T1.txt", "data/Unordered20X20_T1.txt", "data/Ordered20X20_T2.txt", "data/Unordered20X20_T2.txt"]

# nested lists containing values belonging to each file
ϵ = []
m = []

for filename in filenames:
    with open(filename, "r") as file:
        file.readline()
        lines = file.readlines()

        # lists to be filled with values from current file
        ϵ_vals = []
        m_vals = []

        # data is in one column with values separated by ","
        for line in lines:
            values = line.strip().split(",")
            ϵ_vals.append(float(values[0]))
            m_vals.append(float(values[1]))
        
        # adding to nested lists
        ϵ.append(ϵ_vals)
        m.append(m_vals)

steps = len(ϵ[0])
steps = np.arange(steps)

# plot energy data
fig = plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(steps, ϵ[0], 'o', color='orange', label=r'Ordered', markersize=1)
plt.plot(steps, ϵ[1], 'o', color='crimson', label=r'Unordered', markersize=1)
plt.legend(markerscale=5)
plt.xticks([])
plt.ylabel(r'$\left<\varepsilon\right>$ [$\:J$]')
plt.title(r'$T=1.0\:J/k_\text{B}$')
plt.xlim([-500, 10000])
plt.ylim([-2.05, -1.95])

plt.subplot(2, 2, 2)
plt.plot(steps, ϵ[2], 'o', color='orange', label=r'Ordered', markersize=1)
plt.plot(steps, ϵ[3], 'o', color='crimson', label=r'Unordered', markersize=1)
plt.legend(markerscale=5)
plt.xticks([])
plt.title(r'$T=2.4\:J/k_\text{B}$')
plt.xlim([-4500, 90000])
plt.ylim([-1.28, -1.18])

plt.subplot(2, 2, 3)
plt.plot(steps, m[0], 'o', color='skyblue', label=r'Ordered', markersize=1)
plt.plot(steps, m[1], 'o', color='rebeccapurple', label=r'Unordered', markersize=1)
plt.legend(markerscale=5)
plt.ylabel(r'$\left<|m|\right>$')
plt.xlim([-500, 10000])
plt.ylim([0.95, 1.05])

plt.subplot(2, 2, 4)
plt.plot(steps, m[2], 'o', color='skyblue', label=r'Ordered', markersize=1)
plt.plot(steps, m[3], 'o', color='rebeccapurple', label=r'Unordered', markersize=1)
plt.legend(markerscale=5)
plt.xlim([-4500, 90000])
plt.ylim([0.4, 0.5])

fig.supxlabel(r'Number of Monte Carlo cycles')
plt.tight_layout()
plt.show()
plt.show()