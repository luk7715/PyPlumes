# This file is a part of PyPlumes, a Python model for dust ejection dynamics
# on atmosphereless surfaces
# Version 1.0.0
# Translated and modified from the Fortran-95 code in the following paper by Anastasiia Ershova and Jürgen Schmidt, 
# Two-body model for the spatial distribution of dust ejected from an atmosphereless body, 2021, A&A, 650, A186 

# File: distributions.py
# Description: The functions describing the ejection process and auxilary functions used by them

# Author: Eulrika(Yiqi) Wu
# E-mail: ulkw517@g.ucla.com

import numpy as np
import const
import variables as var



def ejection_direction_distribution(distribution_shape, wpsi, psi, lambdaM, zeta, eta):
  '''
  To calculate the axisymmetric distribution of ejection direction
  Args:
    -distribution_shape: Integer. Used to select the expression for the distribution function
    -wpsi: The polar angle in the coordinate system where the distribution is axisymmetrical
    -psi: The polar angle in the horizontal coordinate system
    -lambdaM:The azimuth in the horizontal coordinate system
    -zeta: Zenith angle of the distribution symmetry axis in the horizontal coordinate system
    -eta: Azimuth angle of the distribution symmetry axis in the horizontal coordinate system
  Return:
    -fpsi: Angle between the ejectiong direction and the symmetry axis 
  Return Type:
    Float
  '''

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


    case 2:
    # Uniform distribution of polar angle inside a cone,
    # uniform distribution of azimuth
      if(wpsi <= omega10) :
        fpsi = 1.0 / (1.0 - np.cos(omega10)) / const.twopi
      else:
        fpsi = 0.0
    


    case 3:
    # pseudo Gaussian distribution of polar angle,
    # uniform distribution of azimuth
      if(wpsi < const.halfpi * 0.99) :
        fpsi = np.exp(-(wpsi-psimax45)**2 / 2.0 / omega45 / omega45)
        # this factor is normalization due to the fact that fpsi
        # domain is from 0 to pi/2 and not from -infinity to +infinity
        fpsi = fpsi / normconst3
        fpsi = fpsi / const.twopi
      else:
        fpsi = 0.0



    case 4:
    # HERE IS THE PLACE FOR WRITING YOUR OWN DISTRIBUTION
     fpsi = 0.0

   
    
  if zeta != 0.0 :
    Jpsi = Jacobian_tilt(psi, lambdaM, zeta, eta)
    fpsi = fpsi * Jpsi

  fpsi = fpsi * np.sin(wpsi)


  return fpsi
    

      

def Jacobian_tilt(psi, lambdaM, zeta, A):  
  '''
  To calculate the jacobian of coordinate transformation from the vertical coordinate system to 
  the coordinate system with z-axis coinciding with the jet axis of symmetry.
  Args:
    -psi: Polar angle of ejection in the horizaontal coordinate system
    -lambdaM: Azimuth of ejection in the horizontal coordinate system
    -zeta: Zenith angle of the distribution symmetry axis in the horizontal coordinate system
    -A: Azimuth of the distribution symmetry axis in the horizontal coordinate system
  Return:
    -J: The value of the jacobian
  Return Type:
    Float
  '''

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
  





def ejection_speed_distribution(u, R):
  '''
  This functionrepresents the possibly time-dependent ejection speed distribution.
  When necessary, go to case 3 to create your own distribution functions. 

  Args:
    -u: ejection speed in m/s
    -R: particle radius in microns
  Return:
    -fu:
  Return Type:
    Float
  '''
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
      fu = 0.0   #replace with fu= your desired distribution 
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
  '''
  This function represents the size distribution of the dust particles.
  It can be used also to obtain the mean radius, cross section
  or volume of the dust particles.
  When necessary, go to case 4 to create your own distribution function.

  Args:
   -R: Particle radiusin microns
   -sd: integer number used to select the distribution 
   -mom:   parameter defining the obtained quantity: 
           0 -- number density, 1 -- mean radius, 
           2 -- cross section, 3 -- volume 
  Return:
  -fR
  Return Type:
  -Float 

  '''
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
      q = 5.0
      fR = R**(-3)
      C_size_distr = (r2**(1.0-q) - r1**(1.0-q)) / (1.0 - q)
       #Replace with desired distribution
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
  '''
  This function represents the gamma(t) factor in formula 43 from paper.
  When necessary, go to case 4 to create your own distribution function.

  Args:
   -t: moment of ejection, seconds
   -gamma0:  parameter which can be used in the definition
             of the function. In the implemented examples
             gamma0 is the production rate at maximum
             (particles/second)
   -ratefun: parameter used to choose the expression for the production rate 
  Returns:
   -gammarate
  Return Type:
   Float
  '''
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
      gammarate = 0.0 #replace with desired distribution 
  # endselect
  
  return gammarate
# # end def production_rate

