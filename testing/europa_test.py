import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
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
import distributions as distf
import gu

r1 = 0.20
r2 = 20.0
europa_orbital_period = 3.06822*(10**5)	## in seconds
# number of sources
Ns = 4
#number of points on the SC trajectory for which the number density is to be calculated
nt = 80
tnow = 0.0

#readin the data
dphi, var.point, var.source = input.get_europa_input(Ns,nt)


#integrate the production rate over time to find the total mass
m1 = gu. mass_production(var.source.sd[0], r1, r2)
m2 = gu.mass_production(var.source.sd[1], r1, r2)
mass_steep = var.source.production_rate * m2
mass_shallow = var.source.production_rate * m1

print("with the shallow size distribution, the total mass produced in a second is " + str(mass_shallow) +"kg\n")
print("with the steep size distribution, the total mass produced in a second is " + str(mass_steep) +"kg\n")


#define the parameters for integration
massflux = np.zeros([nt,2])
eadini = var.source.ejection_angle_distr
sdini = var.source.sd
uiini = var.source.ui
GUini = var.source.Gu_precalc
rvector_initial = var.point.rvector 
beta_initial = var.point.beta

#Loop through thetwo sources
for i_s in range(0, Ns):
	
	var.source.ejection_angle_distr = var.source.ejection_angle_distr[i_s]
	var.source.sd = var.source.sd[i_s]
	var.source.ui = var.source.ui[: ,i_s]
	var.source.Gu_precalc = var.source.Gu_precalc[:,i_s] 

#Loop through all data points
	for i in range(0,nt):
		var.point.beta = var.point.beta[i]
		var.point.rvector = var.point.rvector[i,:]

		massflux[i,:] = integrator.DUDI(tnow)

		var.point.rvector = rvector_initial
		var.point.beta = beta_initial


	output.surface_deposition_out(i_s+1, massflux[:,0], nt, dphi)

#reset the parameters 
	var.source.ejection_angle_distr = eadini
	var.source.sd = sdini 
	var.source.ui = uiini 
	var.source.Gu_precalc = GUini
	

	
###plotting the results
fnames = ('PyPlumes/results/narrow_jet_shallow_sd', 'PyPlumes/results/diffuse_source_steep_sd', 'PyPlumes/results/diffuse_source_shallow_sd', 'PyPlumes/results/narrow_jet_steep_sd')
ead = (' narrow jet', ' diffuse source')
sd = ('steep size distribution,', 'shallow size distribution,')
d = []
rad = []
for i in range(4):
	d.append(np.loadtxt(fnames[i]+'.dat', usecols=range(2)))
	rad.append(d[i][:,0])
plt.yscale('log')
plt.xlim(0,110)
plt.plot(rad[0], d[0][:,1], 'k--', label = sd[1] + ead[0])
plt.plot(rad[2], d[2][:,1], 'b--', label = sd[1] + ead[1])
plt.plot(rad[3], d[3][:,1], 'k-', label = sd[0] + ead[0])
plt.plot(rad[1], d[1][:,1], 'b-', label = sd[0] + ead[1])
plt.ylabel('deposited mass ($kg/m^2/s$)')
plt.xlabel('distance from the source (m)')
plt.suptitle('Deposited Mass VS Distance for Different Ejections')
plt.legend(loc = 'upper right')
imname = 'PyPlumes/results/mass_deposition' + '.png'
plt.savefig(imname)
plt.show()
