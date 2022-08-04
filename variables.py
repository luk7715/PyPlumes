# This file is a part of DUDI, the Fortran-95 implementation 
# of the two-body model for dust dynamics
# Version 1.0.0
# This is free software. You can use and redistribute it 
# under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
# If you do, please cite the following paper
# Anastasiia Ershova and Jürgen Schmidt, 
# Two-body model for the spatial distribution of dust ejected from
# an atmosphereless body, 2021, A&A, 650, A186 
# File: define_types.f95
# Description: Definition of the structures used by the routines

# Author: Anastasiia Ershova
# E-mail: vveyzaa@gmail.com

import input

class ejection_speed_properties:     # ejection speed properties 
  def __init__(self, ud_shape, umax, umin):
    self.ud_shape = ud_shape                 # parameter defining which distribution is used for ejection speed
    self.umax = umax                         # gas velocity
    self.umin = umin                         # a parameter in velocity distribution                

#define ud first so it can be used for the following class
ud = ejection_speed_properties()



class source_properties(ejection_speed_properties):                     # parameters of dust ejection
  def __init__(self, rrM, r, alphaM, betaM,zeta,eta,symmetry_axis,ejection_angle_distr,\
    sd,ui,Gu_precalc,production_fun,production_rate,is_jet,ejection_speed_properties):
    self.rrM = rrM                           # Cartesian coordinates of a point source in the moon-centered coordinate system
    self.r = r
    self.alphaM = alphaM                     # polar angle of the point source
    self.betaM = betaM                       # eastern longitude of the point source
    self.zeta = zeta                         # zenith angle of the axis around which ejection is symmetrical
    self.eta = eta                           # azimuth of this axis (counted from the local North, clockwise)
    self.symmetry_axis = symmetry_axis       # unit vector in moon-centered coordinate system pointing to the direction of the axis around which ejection is symmetrical
    
    #type(ejection_speed_properties) ud       # parameters of ejection speed distribution
    self.ud_shape = ejection_speed_properties.ud_shape
    self.ud_umax = ejection_speed_properties.umax
    self.ud_umin = ejection_speed_properties.umin
   
    self.ejection_angle_distr = \
      ejection_angle_distr                   # parameter defining which ejection angle distribution is used; 1 -- Gaussian, 2 -- uniform cone, 3 -- parabola inside a cone
    self.sd = sd                             # parameter to select ejected dust size distribution
    self.ui = ui                             # interpolation grid for GRm(u,Rmin,Rmax) precalculation
    self.Gu_precalc = Gu_precalc             # Gu(Rmin,Rmax)
    self.production_fun = production_fun     # parameter used to select a def for production rate (if <= 0, the production rate is constant)
    self.production_rate = production_rate   # parameter used in definition of the def for production rate (usually the normalization factor)
    self.is_jet = is_jet                     #  True  if the ejection is concentrated (omega is small)
  # end type source_properties

#define the variables!!!
source = source_properties()


class position_in_space:                      # where the dust density is to be calculated
  def __init__(self, r, r_scaled, alpha, beta, rvector, compute):
    self.r = r                                # distance from the center of the moon
    self.r_scaled = r_scaled                  # distance from the center of the moon divided by the moon radius
    self.alpha = alpha                        # polar angle
    self.beta = beta                          # eastern longitude
    self. rvector = rvector                   # Cartesian coordinates of the point 
    self.compute = compute                    # marker that density should or should not be computed in this point                  
  # end type position_in_space


#define the variables!!!
point = position_in_space()
