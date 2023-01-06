# This file is a part of PyPlumes, a Python model for dust ejection dynamics
# on atmosphereless surfaces
# Version 1.0.0
# Translated and modified from the Fortran-95 code in the following paper by Anastasiia Ershova and Jürgen Schmidt, 
# Two-body model for the spatial distribution of dust ejected from an atmosphereless body, 2021, A&A, 650, A186 


# File: const.py
# Description: The fundamental constants and the numerical parameters that are used by other functions

# Author: Eulrika(Yiqi) Wu
# E-mail: ulkw517@g.ucla.com


import numpy as np
import os 


# This file is used to set the fundamental constants and auxilary numbers
# Modify this file to fit for your chosen body!

pi = 3.1415926535897930
halfpi = pi / 2.00
sqrtpi = np.sqrt(pi)
twopi = 2.00 * pi
sqrt2 = np.sqrt(2.00)
deg2rad = pi / 180.00
rad2deg = 180.00 / pi
gravity_constant = 6.674*(10**(-11))   		# m^3 / kg / s^2
  
# parameters defining the moon
moon_mass = 4.8* (10**22) 	# kg
gm = gravity_constant * moon_mass
rm = 3.12*(10**6)  	# meters
vesc = float(np.sqrt((gm * 2.00) / rm))

      	
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
rmin = 0.001
    # upper boundary for def Gu(rmin, rmax), microns
rmax = 2
    
# parameters controling accuracy
    # number of precalculated values of GR(u) integral
    # <=> u/u_gas goes from 0 to 1 with step 1/GRN
GRN = 2000
    # order of integration G(R,u) over R to obtain GR(u)
order_R = 10
    # order of Gaussian quadrature for integration
    # of n(r, alpha, beta, v, theta, lambda) over velocity (v)
    # separately for bound and unbound particles
order_v_el = 10
order_v_hy = 10
    

