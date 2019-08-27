"""
Plots the output from our hydro solver and saves to file.
"""

import numpy as np
import matplotlib.pyplot as plt

import sys

filename = sys.argv[1]

data = np.genfromtxt(filename, delimiter=",").T
positions = np.linspace(0, 1, len(data[0]))

fig, ax = plt.subplots(1, 3, figsize=(8, 3), sharex=True)
ax = ax.flatten()

ax[0].scatter(positions, data[4], s=1.0)
ax[0].set_ylabel("Density")

ax[1].scatter(positions, data[5], s=1.0)
ax[1].set_ylabel("Velocity")
ax[1].set_xlabel("Position")

ax[2].scatter(positions, data[6], s=1.0)
ax[2].set_ylabel("Pressure")

ax[2].set_xlim(0.25, 0.75)

fig.tight_layout()
fig.savefig("output.pdf")
