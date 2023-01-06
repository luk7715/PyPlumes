# This file is a part of PyPlumes, a Python model for dust ejection dynamics
# on atmosphereless surfaces
# Version 1.0.0
# Translated and modified from the Fortran-95 code in the following paper
# Anastasiia Ershova and Jürgen Schmidt, 
# Two-body model for the spatial distribution of dust ejected from
# an atmosphereless body, 2021, A&A, 650, A186 
# File: gu.py
# Description: The subroutines that performe integration over the observed
#              interval of particle sizes

# Author: Eulrika(Yiqi) Wu
# E-mail: ulkw517@g.ucla.com

import numpy as np
import ultilities as ultilities
import const 
import variables as var
import distributions as distf



def Gu_integral(sd):
  '''
  This function calculates the function G^p_u for the array of u-values (stored in ui).
  There are GRN values of ui uniformly distributed between umin and umax

  Args:
   -sd: Integer used to select distribution shape
   '''

  ui = np.zeros(const.GRN)

  for i in range (0,const.GRN):
    ui[i] = var.ud.umin + float(i+1) / float(const.GRN)*(var.ud.umax - var.ud.umin)
  

  up = const.rmax
  low = const.rmin
  Si = np.zeros(const.GRN)
  
  Xi, Wi = ultilities.GaussLegendreQuadra(const.order_R)
  for i  in range(0, const.GRN):
    ldif = up - low
    ldif = ldif * 0.50
    lsum = up + low
    lsum = lsum * 0.50  	

    for ii  in range(0, const.order_R):
      r_i = ldif * Xi[ii] + lsum
      Si[i] = Si[i] + ldif * Wi[ii] \
        * distf.size_distribution(r_i, sd, const.p) \
        * distf.ejection_speed_distribution(ui[i], r_i)
   
  
  # S * ice density [kg/m^3] * 1d-18 m^3/micron^3 * 4/3 pi, this is to obtain mass in kilogramms
  if(const.p == 3):
     Si = Si * const.rho * 1*(10**-18) * 4.0 / 3.0 * const.pi

  # S * pi * 1d-12 [m^2 / micron^2] to obtain area of crossection in m^2
  if(const.p == 2):
     Si = Si * const.pi * 1*(10**-12)
  
  for i  in range(0, const.GRN):
    if(Si[i] < 0  or  Si[i]  !=  Si[i]) :
      f = open("PyPlumes/results/Gu_integral_outputs.txt", "a")
      f.write("\nGu has an incorrect value of " + Si[i] + "\n")
      f.close

  return ui, Si
  



# computes the mass confined between r1 and r2
# resuls is in kg
def mass_production(sd, r1, r2):
  mass = 0.0
  Xi , Wi = ultilities.GaussLegendreQuadra(const.order_R)
  ldif = r2 - r1
  ldif = ldif * 0.50
  lsum = r2 + r1
  lsum = lsum * 0.50  
  
  for i  in range(0, const.order_R):
    r_i = ldif * Xi[i] + lsum
    mass = mass + ldif * Wi[i] * distf.size_distribution(r_i, sd, 3)
  # enddo
  
  # upon integration we obtained the volume of particles as if they were cubes
  # we multiply it by 4pi/3 to obtain the volume of spherical particles
  # 1d-18 is the factor for the conversion from microns^3 to meters^3
  # :, by multiplying by density rho, we obtain the mass
  mass = mass * const.rho * (1*10**-18) * 4.0 / 3.0 * const.pi

  return mass

# # end def mass_production



# end module Gu
