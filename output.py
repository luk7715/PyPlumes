# This file is a part of DUDI, the Fortran-95 implementation 
# of the two-body model for dust dynamics
# Version 1.0.0
# This is free software. You can use and redistribute it 
# under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
# If you do, please cite the following paper
# Anastasiia Ershova and Jürgen Schmidt, 
# Two-body model for the spatial distribution of dust ejected from
# an atmosphereless body, 2021, A&A, 650, A186 
# File: dataoutmod.f95
# Description: The subroutines used to write the result to a text-file

# Author: Anastasiia Ershova
# E-mail: vveyzaa@gmail.com

import const
import variables as var
import numpy as np


def result_out(density,nt):
  
  f = open("PyPlumes/results/twobody_model_result.txt","a") 
  for i  in range(0, nt):
    f.write(str(density[i,:]) + " " + str(density[i,:] )+ " " + str(var.point.r) + \
       " " + str((const.halfpi - var.point.alpha[i]) * const.rad2deg) \
        + " " + str(var.point.beta[i]*const.rad2deg))
    f.write("\n")
  f.close()



# # end def result_out



def cassini_flyby_out(density, ttab, bg, nt):
  #integer, intent(in) :: nt
  #real, intent(in) :: density(nt,2), ttab(nt), bg
  #integer i
  
  f = open("PyPlumes/results/E2_profile.txt","a")
  for i  in range(0, nt):
    densitysum = np.sum (density[i,:]) + bg  
    f.write(str(ttab[i]) + " " + str(densitysum) + "\n")
  f.close 
# # end def cassini_flyby_out  



