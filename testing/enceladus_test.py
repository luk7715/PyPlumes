import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import variables as var
import integrator
import input
import output

# number of sources
Ns = 1
# number of points on the SC trajectory for which the number density is to be calculated
nt = 100
tnow = 0.0
bg = 0.01
fname = "PyPlumes/input_data_files/Enceladus_jet.dat"
density = np.zeros([nt,2])
tmp_res = np.zeros([nt,2])

var.source = input.read_sources_params(fname,Ns)

#loop through all sources in case there are multiple
for i in range (0, Ns):
    for x in range (0,nt): 
        #loop through all data point
        ttab, var.point = input.read_Cassini_E2(nt)      
        var.point.r = var.point.r[x]; 
        var.point.alpha = var.point.alpha[x]
        var.point.beta = var.point.beta[x]; var.point.r_scaled = var.point.r_scaled[x]
        var.point.rvector = var.point.rvector[x,:]
        
        #run the integration functions for each point to get the corresponding density 
        tmp_res[x,:] = integrator.DUDI(tnow)
    density = density + tmp_res

output.cassini_flyby_out(density, ttab, bg, nt)

d = np.loadtxt('Pyplumes/results/E2_profile.dat', usecols=range(2))
t = d[:,0]
dens = d[:,1]
d = np.loadtxt('Pyplumes/input_data_files/E2_1.6.txt', usecols=range(2))
hrdt = d[:,0]
hrddens = d[:,1]
plt.figure(1)
plt.ylabel('number density of grains m$^{-3}$')
plt.xlabel('seconds from the moment of the closest approach between Cassini and Enceladus')
plt.ylim(0,0.12)
plt.suptitle('Number density profile of the plume during Enceladus flyby')
model = plt.plot(t, dens, label = 'model number density', color = "purple")
hrd = plt.scatter(hrdt, hrddens, label = 'HRD measurements')
plt.legend(loc = 'upper left')
plt.show()
