# This file is a part of DUDI, the Fortran-95 implementation 
# of the two-body model for dust dynamics
# Version 1.0.0
# This is free software. You can use and redistribute it 
# under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
# If you do, please cite the following paper
# Anastasiia Ershova and Jürgen Schmidt, 
# Two-body model for the spatial distribution of dust ejected from
# an atmosphereless body, 2021, A&A, 650, A186 
# File: const.f95
# Description: The fundamental constants and the numerical parameters
#              that are used by many subroutines

# Author: Anastasiia Ershova
# E-mail: vveyzaa@gmail.com

##Python version modified by Eulrika Wu 

import numpy as np
import os 


# fundamental constants and auxilary numbers
pi = 3.1415926535897930
halfpi = pi / 2.00
sqrtpi = np.sqrt(pi)
twopi = 2.00 * pi
sqrt2 = np.sqrt(2.00)
deg2rad = pi / 180.00
rad2deg = 180.00 / pi
gravity_constant = 6.674*(10**(-11))   		# m^3 / kg / s^2
  
# parameters defining the moon
moon_mass = 1.08022*(10**(20))  	# kg
gm = gravity_constant * moon_mass
rm = 252*(10**3)  	# meters
vesc = float(np.sqrt((gm * 2.00) / rm))
#print("vesc = " + str(vesc))
      	
# other parameters defining the quantity of interest
    # density of the dust particles (if one wants to compute mass)
    # if the quantity which we want to compute is dust flux through
    # the surface parallel to the moon surface
    # parameter flux should be  True , if it is  False  : density is computed
rho = 920.00  	# kg/m^3 
flux =  False 
    
# parameters of Gu def
    # p = 0 -- number density is computed, 1 -- mean radius,
    # 2 -- cross section, 3 -- mass density
p = 0
    # lower boundary for def Gu(rmin, rmax), microns
rmin = 1.60
    # upper boundary for def Gu(rmin, rmax), microns
rmax = 6.00
    
# parameters controling accuracy
    # number of precalculated values of GR(u) integral
    # <=> u/u_gas goes from 0 to 1 with step 1/GRN
GRN = 2000
    # order of integration G(R,u) over R to obtain GR(u)
order_R = 30
    # order of Gaussian quadrature for integration
    # of n(r, alpha, beta, v, theta, lambda) over velocity (v)
    # separately for bound and unbound particles
order_v_el = 50
order_v_hy = 20
    
#print('\n' + 'Constants and parameters imported!!')

    	
# end module const
