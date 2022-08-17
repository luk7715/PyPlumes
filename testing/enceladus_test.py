import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import const 
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
#print("ui + Si + axis = " + str(var.source.ui)+ str(var.source.Gu_precalc) + str(var.source.symmetry_axis))
#ttab, var.point = input.read_Cassini_E2(nt)

#print (var.source.r)

for i in range (0, Ns):
    for x in range (0,nt): 
        ttab, var.point = input.read_Cassini_E2(nt)
       # print(var.point.rvector)

       
        var.point.r = var.point.r[x]; 
        var.point.alpha = var.point.alpha[x]
        var.point.beta = var.point.beta[x]; var.point.r_scaled = var.point.r_scaled[x]
        var.point.rvector = var.point.rvector[x,:]
       # print(var.point.rvector)

        #var.source = var.source[i]
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
plt.ylabel('number density of grains $> 1.6\\ \\mu m$')
plt.xlabel('seconds from the moment of the closest approach')
plt.ylim(0,0.12)
plt.suptitle('HRD number density profile of E2 flyby')
model = plt.plot(t, dens, label = 'model number density', color = "purple")
hrd = plt.scatter(hrdt, hrddens, label = 'HRD measurements')
plt.legend(loc = 'upper left')
plt.show()
