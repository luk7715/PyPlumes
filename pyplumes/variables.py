##Classes of variables
import numpy as np

class ejection_speed_properties:
  ''' This class defines parameters in ejection speed distribution'''
       # ejection speed properties 
  def __init__(self, ud_shape, umax, umin):
    self.ud_shape = ud_shape                 # integer parameter used to select the desired distribution
    self.umax = umax                         # maximal possible ejection speed, m/s
    self.umin = umin                         # minimal possible ejection speed, m/s

#define ud first so it can be used for the following class
ud = ejection_speed_properties(0,0,0)



class source_properties(ejection_speed_properties): 
  '''
  This class contains parameters describing the dust ejection'''
  def __init__(self, rrM, r, alphaM, betaM,zeta,eta,symmetry_axis,ejection_angle_distr,\
    sd,ui,Gu_precalc,production_fun,production_rate,is_jet,ejection_speed_properties):
    self.rrM = rrM                           # Cartesian coordinates of a point source in the moon-centered coordinate system
    self.r = r
    self.alphaM = alphaM                     # polar angle of the point source
    self.betaM = betaM                       # eastern longitude of the point source
    self.zeta = zeta                         # zenith angle of the axis around which ejection is symmetrical
    self.eta = eta                           # azimuth of this axis (counted from the local North, clockwise)
    self.symmetry_axis = symmetry_axis       # unit vector in moon-centered coordinate system pointing to the direction of the axis around which ejection is symmetrical
    self.ejection_angle_distr = \
      ejection_angle_distr                   # parameter defining which ejection angle distribution is used; 1 -- Gaussian, 2 -- uniform cone, 3 -- parabola inside a cone
    self.sd = sd                             # parameter to select from given ejected dust size distribution
    self.ui = ui                             # interpolation grid for Gu(u,Rmin,Rmax) precalculation, array of GRN elements
    self.Gu_precalc = Gu_precalc             # Precalculated Gu(Rmin,Rmax), array of GRN elements
    self.production_fun = production_fun     # parameter used to select a def for production rate (if <= 0, the production rate is constant)
    self.production_rate = production_rate   # parameter used in definition of the def for production rate (usually the normalization factor)
    self.is_jet = is_jet                     #  True  if the ejection is concentrated (recommendo the define as True when source.omega < 0.1)
   
    #type(ejection_speed_properties) ud       # parameters of ejection speed distribution
    self.ud_shape = ejection_speed_properties.ud_shape
    self.ud_umax = ejection_speed_properties.umax
    self.ud_umin = ejection_speed_properties.umin
   
   
  # end type source_properties


source = source_properties(0,0,0,0,0,0,0,0,0,0,0,0,0,0,ud)



class position_in_space:  
  '''
  This class defines the point in space where the density is calculated
  Measured in the moon's centered coordinate system
  '''
  def __init__(self, r, r_scaled, alpha, beta, rvector,compute):
    self.r = r                                # Radial distance from the center of the moon to the spacecraft, meters
    self.r_scaled = r_scaled                  # Radial distance from the center of the moon to the spacecraft in units of the moon radius
    self.alpha = alpha                        # Polar angle of the spacecraft, in radians 
    self.beta = beta                          # Eastern longitude of the spacecraft, in radians
    self.rvector = rvector                    # Cartesian coordinates of the spacecraft
    self.compute = compute                    # TRUE if the density in this point must be calculated (allows to exclude regions where the calculations are unnecessary)
                                                          
  # end type position_in_space

point = position_in_space

pointS = np.array([[ [position_in_space for col in range(128)] for col in range(128)] for row in range(128)])


