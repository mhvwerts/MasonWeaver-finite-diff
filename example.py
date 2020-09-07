#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example of using the finite-difference Mason-Weaver solver.

M. H. V. Werts, 2020
CNRS/ENS Rennes, lab. SATIE

published under the CeCILL v2.1 license
    https://cecill.info/licences/Licence_CeCILL_V2.1-en.txt
"""

import numpy as np
import matplotlib.pyplot as plt

from masonweaver_FDsolver import MWsoln


Ntimesteps = 500
Ngridpoints = 1000


# Calculate the finite-difference solution of the Mason-Weaver equation
# on a grid of 'Ngridpoints' grid points
# over 'Ntimesteps' time steps (exponentially increasing delta t)
# The parameters 'D' and 'sg' are for 60 nm diam. spherical gold nanoparticles
# in water at +4°C.

soln = MWsoln(z_max = 0.01, t_end = 3.6E6,
              J = Ngridpoints, N = Ntimesteps,
              D = 4.32E-12, 
              sg = 2.47E-9*9.81) # for 60 nm gold nanoparticle at 4C


# Plot/analyze the solution
# (This demonstrates how to use the 'MWsoln' object containing the solution)

plt.figure(1)
plt.clf()
plt.title('Sedimentation 60nm diam. gold NP in water at 4°C')
for i in range(0, Ntimesteps, 50):
    plt.plot(soln.z[:-1]*100, soln.c[:-1,i],
             label='t={:.0f} h '.format(soln.t[i]/3600))
plt.ylim(ymin = 0.0, ymax = 3.0)
plt.xlabel('height / cm')
plt.ylabel('rel. concentration')
plt.legend()

# time behaviour
plt.figure(2)
plt.clf()
plt.title('Concentration at bottom of cell')
plt.plot(soln.t/3600, soln.c[1,:].T)
plt.ylabel('rel. concentration')
plt.xlabel('time / h')

# convergence to analytical steady-state
plt.figure(3)
plt.clf()
plt.title('Sedimentation equilibrium gradient')
B = soln.z_max / (soln.z0*(1-np.exp(-soln.z_max/soln.z0)))
plt.plot(soln.z*100, soln.c[:,-1],
         label = 'finite diff.')
plt.plot(soln.z*100, B * np.exp(-soln.z/soln.z0),
         label = 'analytic')
plt.xlabel('height / cm')
plt.ylabel('rel. concentration')
plt.legend()

plt.show()

