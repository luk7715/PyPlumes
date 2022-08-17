# This file is a part of DUDI, the Fortran-95 implementation 
# of the two-body model for dust dynamics
# Version 1.0.0
# This is free software. You can use and redistribute it 
# under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
# If you do, please cite the following paper
# Anastasiia Ershova and Jürgen Schmidt, 
# Two-body model for the spatial distribution of dust ejected from
# an atmosphereless body, 2021, A&A, 650, A186 
# File: distributions_fun.f95
# Description: The functions describing the ejection process and auxilary
#              functions used by them

# Author: Anastasiia Ershova
# E-mail: vveyzaa@gmail.com

import numpy as np
import const
import variables as var


# The axisymmetric distribution of ejection direction
# distribution_shape is the parameter used to select
# the expression for the PDF
# wpsi is the polar angle in the coordinate system where
# the distribution is axisymmetrical
# psi is the polar angle in the horizontal coordinate system
# lambdaM is the azimuth in the horizontal CS
# zeta and eta are respectively zenith angle and azimuth
# of the distribution symmetry axis in the horizontal CS
def ejection_direction_distribution(distribution_shape, wpsi, psi, lambdaM, zeta, eta):

###Define constants 
  normconst1 = 7.5960829056967811*(10**-3)
  normconst3 = 8.5760756217641998*(10**-1)
  psimax0 = 0.0
  psimax45 = 0.7853982
  omega3 = 0.05235988
  omega5 = 0.08726646
  omega10 = 0.1745329
  omega45 = 0.7853982
  fpsi = float(0.0)
  Jpsi = float()

  match distribution_shape:
    case 1:
    # pseudo Gaussian distribution of polar angle, uniform distribution of azimuth
      if(psi < const.halfpi * 0.99  and  wpsi < const.halfpi * 0.99) :
        fpsi = np.exp((-(wpsi-psimax0)**2) / 2.0 / omega5 / omega5)
        # this factor is normalization due to the fact that fpsi
        # domain is from 0 to pi/2 and not from -infinity to +infinity
        fpsi = fpsi / normconst1
        fpsi = fpsi / const.twopi
      else:
        fpsi = 0.0
      # endif

    case 2:
    # Uniform distribution of polar angle inside a cone,
    # uniform distribution of azimuth
      if(wpsi <= omega10) :
        fpsi = 1.0 / (1.0 - np.cos(omega10)) / const.twopi
      else:
        fpsi = 0.0
      # endif


    case 3:
    # pseudo Gaussian distribution of polar angle,
    # uniform distribution of azimuth
      if(wpsi < np.halfpi * 0.99) :
        fpsi = np.exp(-(wpsi-psimax45)**2 / 2.0 / omega45 / omega45)
        # this factor is normalization due to the fact that fpsi
        # domain is from 0 to pi/2 and not from -infinity to +infinity
        fpsi = fpsi / normconst3
        fpsi = fpsi / const.twopi
      else:
        fpsi = 0.0
      # endif


    case 4:
    # HERE IS THE PLACE FOR WRITING YOUR OWN PDF
     fpsi = 0.0

   
    
  if zeta != 0.0 :
    Jpsi = Jacobian_tilt(psi, lambdaM, zeta, eta)
    fpsi = fpsi * Jpsi
  # endif

  fpsi = fpsi * np.sin(wpsi)

  #print("fpsi = "+ str(fpsi))

  return fpsi
    
# # end def ejection_direction_distribution


      

# Jacobian of coordinate transformation from vertical CS
# to the CS with z-axis coinciding with the jet axis of symmetry
# zeta and A are respectively zenith angle ang azimuth
# of the distribution symmetry axis in the horizontal CS
# psi and lambdaM are respectively polar angle
# and azimuth of ejection in the horizontal CS
def Jacobian_tilt(psi, lambdaM, zeta, A):  
  J= 0.0
  if(psi < 0.120  and  zeta < 0.120) :
    J = psi / np.sqrt(psi*psi + zeta*zeta - 2.0 * psi * zeta * np.cos(lambdaM - A))
  else:
    sinpsi = np.sin(psi) ; cospsi = np.cos(psi)
    sinzeta = np.sin(zeta) ; coszeta = np.cos(zeta)

    J = 4.0 * sinpsi / np.sqrt(10.0 - 2.0 * (cospsi - sinpsi) * (cospsi + sinpsi) \
      - 3.0 * np.cos(2.0 * (psi - zeta)) - 2.0 * (coszeta - sinzeta) * (coszeta + sinzeta) \
      - 3.0 * np.cos(2.0 * (psi + zeta)) \
      - 8.0 * np.cos(2.0 * (lambdaM - A)) * sinpsi * sinpsi * sinzeta * sinzeta \
      - 32.0 * np.cos(lambdaM - A) * sinpsi * cospsi * sinzeta * coszeta)
  # endif

  return J
  
# # end def Jacobian_tilt



# This def represents the ejection speed distribution
# (possibly time-dep# endent)
# ud is  the parameter used to select the expression for the distribution
# u is the ejection speed
# R is the particle size
def ejection_speed_distribution(u, R):
  Rc = 0.5
  fu = 0.0

  match (var.ud.ud_shape):

    case 1 :
      Rrel = R / Rc
      urel = u / var.ud.umax
      fu = Rrel * (1.0 + Rrel) *( (1.0 - urel)**(Rrel - 1.0) )* urel / var.ud.umax
    case 2:
      fu = 0.0
      if(u < var.ud.umax  and  u > var.ud.umin):
         fu = 1.0 / (var.ud.umax - var.ud.umin)

    case 3:
    # HERE IS THE PLACE FOR WRITING YOUR OWN PDF
      fu = 0.0
  # endselect

  return fu
# # end def ejection_speed_distribution




# This def represents the size distribution of the dust particles.
# It can be used also to obtain the mean radius, cross section
# or volume of the dust particles.
# R is a particle radius
# sd is a parameter used to select the expression for the distribution
# mom defines the obtained quantity: 0 -- number density,
# 1 -- mean radius, 2 -- cross section, 3 -- volume
def size_distribution(R, sd, mom):
  mu = -1.0
  sigma = 1.50
  r1 = 0.20
  r2 = 20.0
  fR = 0.0 
  match sd:

    case 1:
      fR = np.exp(-(np.log(R) - mu)**2 / 2.0 / (sigma**2)) / R
      C_size_distr = sigma * const.sqrtpi * np.sqrt(2.0)

    case 2:
      q = 3.0
      if(r1 <= R  and  R <= r2) :
        fR = R**(-3)
      else:
        fR = 0.0
      # endif
      C_size_distr = (r2**(1.0-q) - r1**(1.0-q)) / (1.0 - q)

    case 3:
      q = 5.0
      if(r1 <= R  and  R <= r2) :
        fR = R**(-3)
      else:
        fR = 0.0
      # endif
      C_size_distr = (r2**(1.0-q) - r1**(1.0-q)) / (1.0 - q)

    case 4:
      # HERE IS THE PLACE FOR WRITING YOUR OWN PDF
      fR = 0.0
  # endselect
  fR = ( R**mom )* fR / C_size_distr

  return fR
 # # end def size_distribution




# This def represents the factor gamma(t) in Formula 43
# t is the moment of ejection
# gamma0 is a parameter which can be used in the definition of
# the def. In the implementet examples gamma0 is the production
# rate at maximum.
# ratefun is the parameter used to choose the expression for the gammarate
def production_rate(t, gamma0, ratefun): 
  gammarate = 0.0
  match (ratefun):
    case 1:
      gammarate = 0.0
      if(t < tmax and  t >= 0.0):
        gammarate = gamma0

    case 2:
      gammarate = 0.0
      tmax = 500.00
      if(t > 0.0  and  t < 2.0 * tmax) :
        gammarate = gamma0 * (-t**2 + 2.0 * t * tmax) / (tmax**2)
      # endif
    case 3:
      # HERE IS THE PLACE TO WRITE YOUR OWN def FOR THE PRODUCTION RATE
      gammarate = 0.0
  # endselect
  
  return gammarate
# # end def production_rate

