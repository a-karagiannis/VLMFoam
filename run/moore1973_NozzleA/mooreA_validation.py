#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 08:00:27 2019

@author: akaragiannis
"""

import numpy as np
import matplotlib.pyplot as plt

filename_exp = "mooreA_centerline_exp.csv"
filename_sim = "mooreA_centerline_sim.csv" 


# Set normalisation parameters
P_0 = 25000.0 # [Pa], stagnation pressure at inlet
x_length = 0.75 # [m], nozzle length

# Read experimental and simulation data from respective files
data_exp = np.transpose(np.genfromtxt(filename_exp, delimiter=','))
data_sim = np.transpose(np.genfromtxt(filename_sim, 
                                      skip_header = 1, delimiter=','))

# Assign experimental pressure data and x coordinate to separate arrays
experiment_P_array = data_exp[1, :]
experiment_x_array = data_exp[0, :]

# Same for simulation data
# Adjust simulation data to match format of the Moore et al (1973) paper plots
            # -0.25 to correct for different coordinate starting points
simulation_P_array = data_sim[0, :]/P_0
simulation_x_array = (data_sim[4, :] - 0.25)/x_length 

# Assign the radius data from the simulation
simulation_radius_array = data_sim[1, :]


###############
### PLOTTING
###############

fig, ax1 = plt.subplots()

# Color of data in the plot
colors = ['k'] # Black

# Scatter plot - experimental pressure
scatter_pressure = ax1.scatter(experiment_x_array, experiment_P_array, 
                               marker='s', c=colors, s=40, alpha=1.0, 
                               label = 'Experiment ($P/P_{0}$)')

# Line plot - simulated pressure
plot_pressure, = ax1.plot(simulation_x_array, simulation_P_array, 
                         color = 'k', label = 'Present model ($P/P_{0}$)')

# Plot frame
ax1.set_xlim((-0.25, 0.6))
ax1.set_ylim((0.1, 1.0))
x0,x1 = ax1.get_xlim()
y0,y1 = ax1.get_ylim()

# Text box on figure's upper left corner
props = dict(boxstyle='round', facecolor='white', alpha=1.0)
textstr = '\n'.join((
    r'Moore Nozzle Type A',
    r'$P_{0} = 25 $ kPa',
    r'$T_{0} = 354.6 $ K',
    r'Supersonic outlet'))
ax1.text(0.03, 0.97, textstr, transform=ax1.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

# Axes labels
ax1.set_xlabel('$X/X_{Length}$ [-]', fontsize = 14)
ax1.set_ylabel('$P/P_{0}$ [-]', fontsize = 14)
ax1.tick_params(axis = 'x', labelsize = 12, which = 'both', 
                width = 1.0, length = 5.0)
ax1.tick_params(axis = 'y', labelsize = 12, which = 'both', 
                width = 1.0, length = 5.0)


# Grid 
ax1.grid(b=True, which='major', color='k', linestyle='--', alpha = 0.2)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
plot_radius, = ax2.plot(simulation_x_array, simulation_radius_array, 
                        color = 'k', linestyle="--", 
                        label = 'Present model ($r_{d}$)')
scatter_radius = ax2.scatter(0.49308, 0.025e-06, marker='^', c=colors, 
                             s=40, alpha=1.0, label = 'Experiment ($r_{d}$)') 
                                                    # Measured droplet radius
ax1.set_adjustable("datalim")

# Switch to ax2 for the radius plot on the right vertical axis 
ax2.set_yscale('log')  # Logarithmic scale
ax2.set_xlim((-0.25, 0.6))
ax2.set_ylim((1e-10, 1e-6))
x2,x3 = ax2.get_xlim()
y2,y3 = ax2.get_ylim()

ax2.set_ylabel('Droplet radius [m]', fontsize = 14)
ax2.tick_params(axis = 'y', labelsize = 12, which = 'major', 
                width = 1.0, length = 5.0)

# Legend on figure's upper right corner
plots = [plot_pressure, plot_radius, scatter_pressure, scatter_radius]
ax1.legend(plots, [plots_.get_label() for plots_ in plots], 
           loc= 'upper right', edgecolor='k', fontsize = 10)

fig.tight_layout()  # otherwise the right y-label is slightly clipped

plt.show()

# Output figure as scaled .pdf
plt.savefig('Validation.pdf', bbox_inches = 'tight', pad_inches = 0)  
