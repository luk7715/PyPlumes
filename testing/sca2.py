import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys
import inspect
from scipy import integrate as integrate

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import variables as var
import input
import output
import distributions as distf
import const
import los
import integrator as inte

r1 = 0.0010 #minimum dust size in micron
r2 = 2.0 #maximum dust size in micron
rp = 1000

mass_rate_ky = 1.5 #mass production rate from Keat Yeoh in kg/s
dustsize= np.array([10.0, 50.0, 100.0, 500.0, 1000.0])*(10**-3) #dust size values from Keat Yeoh et al poster, micron
production_rate_ky = mass_rate_ky/const.rho /((dustsize**3 )* (10**-18)) # dust production rate in number of particles/s
print(production_rate_ky)

europa_orbital_period = 3.06822*(10**5)	## in seconds
Ns = 1 # number of sources
nt = 5 #number of points on the SC trajectory for which the number density is to be calculated
tnow = 0.0

#readin the data
dphi,var.point,var.source = input.europascatter_input2(Ns,nt)

number = 6.282*(10**12)*3.12*(10**6)*3.12*(10**6)* (5.42005352*(10**-11)+3.99251925*(10**-9))

density = 0.0
tmp_res = 0.0	
# form the integration grid for this specific plume
import numpy as np
from scipy import integrate as integrate

### From the definition of the LOS integral
def LOS_integration(rs,r_vir,r_p):  #### radius in kpc
    rho = lambda l: np.exp(1 - np.log(1+np.sqrt(l**2 + r_p**2)/rs)/(np.sqrt(l**2 + r_p**2)/rs))
    result = integrate.quad(rho,0,np.sqrt(r_vir**2-r_p**2),epsabs=1.49e-08, epsrel=1.49e-08)
    return result[0]

integration_vec = np.vectorize(LOS_integration)  


Vrad = 4.5* 10** (-2)
Hrad = 4.5* 10** (-2)
Vpix = 64
Hpix = 64
ni = 20
sampdist_large = 3**2 #meters
numbertot = 0.0

for i in range(-Hpix, Hpix):
	for ii in range(-Vpix, Vpix):
#computing the 2nd moment of the dust density
#in the nodes along the line of sight
		for iii in range(-ni, ni):		#if the point isn't on the LOS crossing the moon and
		#not on the LOS having 0-index like (i,0) or (0,iii)
			if(var.pointS[i,ii,iii].compute == True):
				var.point = var.pointS[i,ii,iii]
				numbernow = inte.DUDI(tnow)
				print(str('numbernow is ' + str(numbernow)))
			else:
				numbernow = 0.0
			numbertot = numbertot + numbernow

print(str("the 1st total number den is ") + str(numbertot) )
