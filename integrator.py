# This file is a part of DUDI, the Fortran-95 implementation 
# of the two-body model for dust dynamics
# Version 1.0.0
# This is free software. You can use and redistribute it 
# under the terms of the GNU General Public License 
# (http://www.gnu.org/licenses/)
# If you do, please cite the following paper
# Anastasiia Ershova and Jürgen Schmidt, 
# Two-body model for the spatial distribution of dust ejected from
# an atmosphereless body, 2021, A&A, 650, A186 
# File: integrator.f95
# Description: The subroutines that manage the numerical integration

# Author: Anastasiia Ershova
# E-mail: vveyzaa@gmail.com

import numpy as np
import const
import variables as var
import help
import twobody as tb

# change of veps is not recommended
# difference of the actual minimal velocity value
# from the left boundary used in integration
veps = 1*(10**-10)
trpower = 4


# The def rules the integration determining the parameters and the limits
def DUDI(tnow):
  #integer Nprestep
  # fraction which the interval of the pole integration constitutes
  # to the total interval of integration
  prestep_relative_size = 1*(10**-4)
  #integer i
  #real, intent(out) :: density(2)
  #real(8), intent(in) :: tnow
  #real(8) vc, vmax, v_limits(2), amin, vmin, vmin2, amin2
  #type(position_in_space), intent(in) :: point
  #type(source_properties), intent(in) :: source
  #real(8) ldif, lsum, term, vi, f1, f2, vinterval, viprev
  #real(8) angle, dphi, dbeta
  #real xel(order_v_el), wel(order_v_el)
  #real xhy(order_v_hy), why(order_v_hy)
  #logical pole
  
  dphi, dbeta, angle = help.ApuTrajectory() 
  v_limits = np.array([0.0,0.0])
  Nprestep = int(0)

  
  # escape velosity at distance rr
  vc = np.sqrt(20 * const.gm / var.point.r)
  
  amin = 0.0
  amin = (var.point.r + var.source.r) / 4.0 \
      + 0.50 * np.sqrt((var.point.r**2 + var.source.r**2) \
      / 4.0 - var.point.r * var.source.r * np.cos(dphi) / 2.0)

  amin2 = 0.0
  amin2 = (2.0 / var.source.r - var.source.ud_umin**2 / const.gm)**(-1)
  
  if(amin2  !=  amin2):
     amin2 = 0.0

  # minimal velocity possible at given position (defines minimal energy
  # or "size" of the orbit on which a particle can go from rm to rr
  vmin = np.sqrt(np.abs((const.gm * (2.0 / var.point.r - 1.0 / amin))))
  vmin2 = np.sqrt(np.abs(const.gm * (2.0 / var.point.r - 1.0 / amin2)))
  #print("vmin sqrt is " + str(const.gm * (2.0 / var.point.r - 1.0 / amin2)))
  # v_max is : the maximal *possible* speed at radius r,
  # assuming that the ejection velocity is limited by gas velocity

  vmax = np.sqrt(var.source.ud_umax * var.source.ud_umax \
          + 2.0 * const.gm * (1.0 / var.point.r - 1.0 / var.source.r))

  #print("amin is " +str(amin))
  #print(amin2)

  if amin > amin2:
    pole = 1 

  if pole != True:
     vmin = vmin2
  
  density = np.array([0.0,0.0])
  # if the maximal possible velocity is enough to get from rm to rr
  if(vmax > vmin) :
    
    if pole == True :
      Nprestep = estimate_N_steps_pole_integration(var.source.zeta * const.rad2deg, var.point.r_scaled, \
            angle, var.source.is_jet)
      
      #print(Nprestep)

      v_limits[0] = vmin + veps
      v_limits[1] = (vmax - vmin) * prestep_relative_size + vmin
      vinterval = (v_limits[1] - v_limits[0])
      
      f1 =  tb.Integrand_number_density(v_limits[0], amin, dphi, dbeta, tnow)
      viprev = v_limits[0]

      #print(Nprestep)

      for i  in range(1, Nprestep):
        vi = (float(i-1) / float(Nprestep))**trpower * vinterval + v_limits[1]
        f2 = tb.Integrand_number_density(vi, amin, dphi, dbeta, tnow)
        density[0] = density[1] + (vi - viprev) * 0.50 * (f1 + f2)
        f1 = f2
        viprev = vi
      # enddo
    else:
      v_limits[1] = vmin
    # endif
    
    if(vc > v_limits[1]) :  		
      v_limits[0] = v_limits[1]
      v_limits[1] = min(vc, vmax)
      
    # the particles on the elliptic orbits 
      xel,wel = help.GaussLegendreQuadra(const.order_v_el)
      ldif = v_limits[1] - v_limits[0]
      ldif = ldif * 0.50
      lsum = v_limits[1] + v_limits[0]
      lsum = lsum * 0.50
      for i  in range(0, const.order_v_el):
        term = tb.Integrand_number_density(ldif * xel[i] + lsum, amin, dphi, dbeta, tnow)
        density[0] = density[0] + ldif * wel[i]* term
      # enddo
    # endif
      
    if(vc < vmax) :
    # the particles on the escaping trajectories
      xhy ,why = help.GaussLegendreQuadra(const.order_v_hy)
      v_limits[0] = v_limits[1]
      v_limits[1] = vmax
      ldif = v_limits[1] - v_limits[0]
      ldif = ldif * 0.50
      lsum = v_limits[1] + v_limits[0]
      lsum = lsum * 0.50
      for i  in range(0, const.order_v_hy):
        term = tb.Integrand_number_density(ldif * xhy[i] + lsum, amin, dphi, dbeta, tnow)
        density[1] = density[1] + ldif * why[i] * term
      # enddo
    # endif
    # factor indep# endent on velocity
    density = density / var.point.r / var.source.r / np.sin(dphi)
  # endif
  return density
# # end def DUDI

    


def estimate_N_steps_pole_integration(z, r, xi, isjet):

  ximin = 0.17453290  	# 10 degree in radians
  ximax = 0.78539820  	# 45 degree in radians
  Nprestep = 15


  if(isjet == 1) :
    if(ximin < xi, xi < ximax) :
      if(r < 2.0) :
        Nprestep = 15 + 10 * int(z)
        
      else:
        Nprestep = 10 + 5 * int(z)
        
      # endif
    else:
      if(r < 1.05) :
        Nprestep = 80
        
      else:
        Nprestep = 0
        
      # endif
    # endif
 
  else:
    Nprestep = 15
  # endif
  
  #print("n is " + str(Nprestep))
  return Nprestep

# # end def estimate_N_steps_pole_integration




# end module integrator

