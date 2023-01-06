# This file is a part of PyPlumes, a Python model for dust ejection dynamics
# on atmosphereless surfaces
# Version 1.0.0
# Translated and modified from the Fortran-95 code in the following paper
# Anastasiia Ershova and Jürgen Schmidt, 
# Two-body model for the spatial distribution of dust ejected from
# an atmosphereless body, 2021, A&A, 650, A186 
# File: twobody.py
# Description: The functions used to compute the integrand

# Author: Eulrika(Yiqi) Wu
# E-mail: ulkw517@g.ucla.com

import numpy as np
import ultilities as ultilities 
import const
import variables as var
import distributions as distf

##used kwargs or redefine all the parameters 


def Apu_u_angles_ddphidtheta(v,theta,e,dbeta,dphi,u,dphi_is_large): 
  ##access the variables 
  point_alpha = var.point.alpha; point_r = var.point.r; point_beta = var.point.beta
  point_r_scaled = var.point.r_scaled 
  source_alphaM = var.source.alphaM; source_r = var.source.r; source_betaM = var.source.betaM
  source_zeta = var.source.zeta; source_eta = var.source.eta 

  psi = float(); wpsi = float()
  
  #print(point_alpha)

  ##trig calculations 
  sinal = np.sin(point_alpha) ; cosal = np.cos(point_alpha)
  sinalM = np.sin(source_alphaM) ; cosalM = np.cos(source_alphaM)
  sindphi = np.sin(dphi) ; cosdphi = np.cos(dphi)
  sintheta = np.sin(theta)
  
  #	 angular momentum (eq 27)		
  hh = point_r * v * sintheta  	
    #  psi
  psi = np.arcsin(hh / source_r / u)
  #print(psi)
  
  #########################################
  if psi != psi or (np.abs(psi-const.halfpi) < 1*(10**-8)):

    if psi != psi:
      f = open("PyPlumes/results/Apu_angles.txt", "a")
      f.write("\nApu << sin(psi) = " + str(hh/source_r/u)+ " corrections applied\n")
      f.close()
      # endif
    if np.abs(psi-const.halfpi) < 1*(10**-8):
      f = open("PyPlumes/results/Apu_angles.txt", "a")
      f.write("\npsi is close to pi/2, corrections applied\n")
      f.close()
      # endif
    psi = const.halfpi - 1*(10**-5)
  # endif
  ##########################################
  
  # lambda and lambdaM are to be found from a spherical triangle
  # the direction and the length of arc along which a particle
  # traveled from rM to r matters for exact geometry of
  # the spherical triangle namely, sings of lambda and lambdaM
  # dep# end on it
  if source_betaM > point_beta :
    if (source_betaM - point_beta) < const.pi :
      sinlambda = -sinalM * np.sin(dbeta) / sindphi
      coslambda = (cosal * cosdphi - cosalM) / sinal / sindphi
    else:
      sinlambda = sinalM * np.sin(dbeta) / sindphi
      coslambda = (cosal * cosdphi - cosalM) / sinal / sindphi
    # endif  			
  else:
    if(point_beta - source_betaM) < const.pi :
      sinlambda = sinalM * np.sin(dbeta) / sindphi
      coslambda = (cosal * cosdphi - cosalM) / sinal / sindphi
    else:
      sinlambda = -sinalM * np.sin(dbeta) / sindphi
      coslambda = (cosal * cosdphi - cosalM) / sinal / sindphi
    # endif
  # endif

  #print(sinal)
  #print(cosal)

  sinlambdaM = sinal * sinlambda / sinalM
  coslambdaM = (cosal - cosalM * cosdphi) / sinalM / sindphi

  
  if dphi_is_large == True :
    sinlambda = - sinlambda
    coslambda = - coslambda
    sinlambdaM = - sinlambdaM
    coslambdaM = - coslambdaM
  # endif

  var_lambda = float()
  lambdaM = float()

  lambdaM = ultilities.myatan1(coslambdaM, sinlambdaM)
  var_lambda = ultilities.myatan1(coslambda, sinlambda)

  #print("lm = " + str(lambdaM))
  wpsi = np.arccos(np.cos(psi) * np.cos(source_zeta)+ \
    np.cos(lambdaM - source_eta) * np.sin(psi) * np.sin(source_zeta))

  wrr = point_r_scaled ; wvv = float(v / const.vesc)

  #print("v is " + str(v) +"and wvv is " + str(wvv))

  pp = 2.0 * wrr * wrr * wvv * wvv * (sintheta**2)
  dpp = 2.0 * pp / np.tan(theta)
  dee = (wvv * wvv - 1.0 / wrr) * dpp / e
  
  cosphim = (pp - 1.0) / e
  cosphi = (pp / wrr - 1.0) / e

  ddphidtheta = (((1.0 + pp) * (wrr * wvv * wvv - 1.0) + wrr) * 2.0 * np.cos(theta)) \
    / (wvv * np.sqrt(wrr * (-1.0 + wrr + wrr * (wvv**2) - (wrr**3 )* (wvv**2) * (sintheta**2)))) \
    - (2.0 * (sintheta**2 )* (-2.0 + 2.0 * wrr * (wvv**2) + 1.0 / (sintheta**2)))

  ddphidtheta = ddphidtheta * wvv*wvv * wrr / e/e

  numder = float(0)

  if(ddphidtheta  !=  ddphidtheta) :
    delta = 1*(10**-3)
    dphi1 = deltaphi((theta + 2.0*delta), wrr, wvv)
    dphi2 = deltaphi((theta + delta), wrr, wvv)
    dphi3 = deltaphi((theta - delta), wrr, wvv)
    dphi4 = deltaphi((theta - 2.0*delta), wrr, wvv)
    
    numder = (-dphi1 + 8.0 * dphi2 - 8.0 * dphi3 + dphi4) / 12.0 / delta 
    numder = (dphi2 - dphi3) / 2.0 / delta

    f = open("PyPlumes/results/Apu_angles.txt", "a")
    f.write("\nthe derivative d_Delta_phi/d_theta was obtained numerically because the analytical expression contains numerically difficult parts\n")
    f.close()
    ddphidtheta = numder
  # endif  		
  #print(ddphidtheta)
  #print("numder = " + str(numder))
  return psi, wpsi, lambdaM, var_lambda, ddphidtheta, sindphi

  
# # end def Apu_u_angles_ddphidtheta



# compute \Delta\phi from \theta, velocity and spacecraft position
def deltaphi(theta, wrr, wvv):

  deltaphi = float()
  
  pp = 2.0 * (wrr**2) * (wvv**2)* (np.sin(theta)**2)
  e = np.sqrt(1.0 + 2.0 * pp * (wvv * wvv - 1.0 / wrr))
  
  cosphim = (pp - 1.0) / e
  if np.abs(cosphim) > 1.0:
    cosphim = 1.0 * (np.abs(cosphim)/cosphim)


  cosphi = (pp / wrr - 1.0) / e
  if np.abs(cosphi) > 1.0:
    cosphi = 1.0 * (np.abs(cosphi)/cosphi)


  if(theta < const.halfpi) :
    deltaphi = np.arccos(cosphi) - np.arccos(cosphim)
  else:
    deltaphi = 2.0 * const.pi - np.arccos(cosphi) - np.arccos(cosphim)
  # endif
  return deltaphi

# # end def deltaphi



# Integrand_number_density performs integration over theta and lambda,
# returns the expression standing under the integral over v
def Integrand_number_density(velocity, amin, dphi, dbeta, tnow):
 
  semi_major_axis = (2.0 / var.point.r - (velocity**2) / const.gm)**(-1)
  theta = np.array([-999.0,-999.0])
  ee = np.array([0.0,0.0])
  deltat = np.array([0.0,0.0])
  dphi_is_large =  np.array([0,0], dtype = bool) 
  Integrand = float()

  if var.source.production_fun > 0: 
    timeDependence = True
  else: timeDependence = False
  

  # find solutions for theta

  if(semi_major_axis < 0.0  or  semi_major_axis >= amin) :
    if(semi_major_axis > 0.0) :
      ee, theta, deltat, dphi_is_large = theta_geometry_ellipse(var.point.r, var.source.r, velocity,\
         dphi, semi_major_axis, timeDependence)
    else:
      ee, theta, deltat, dphi_is_large = theta_geometry_hyperbola(var.point.r, var.source.r, velocity,\
         dphi, np.abs(semi_major_axis), timeDependence)
 # print(semi_major_axis < 0.0 )
 # print(ee)
#~ 			write(*,*) 'intergand << theta', theta

  uu = np.sqrt(const.vesc* const.vesc + 2.0 * (velocity * velocity / 2.0 - (const.gm / var.point.r)))
  
  if(uu > var.source.ud_umax  or  uu < var.source.ud_umin) :
    fac1 = 0.0
  else:
#	 interpolate Gu from a precalculated table
    #print("ui is " +str(len(var.source.ui)))
    #print(var.source.Gu_precalc)
    if(uu / var.source.ud_umax > var.source.ui[const.GRN-1]) :
      fac1 = velocity * var.source.Gu_precalc[const.GRN-1] / uu / uu
    else:
      fac1 = ultilities.LiNTERPOL(const.GRN, var.source.Gu_precalc, var.source.ui, uu) 
      fac1 = velocity * fac1 / uu / uu
    # endif
  # endif
  
  ####declare variables 
  Integrand = 0.0
  wpsi = np.array([const.halfpi, const.halfpi])
  rate = np.array([0.0,0.0])
  psi = np.array([0.0,0.0])
  ddphidtheta = np.array([0.0,0.0])

  for i  in range(0, 2):
    if(timeDependence == True) :
      rate[i] = distf.production_rate(tnow - deltat[i], var.source.production_rate, \
                                                    var.source.production_fun)
    else:
      rate[i] = var.source.production_rate
    # endif

    if theta[i] >= 0.0  and  rate[i] > 0 :
      
      psi[i], wpsi[i], lambdaM, var_lambda, ddphidtheta[i], sindphi = Apu_u_angles_ddphidtheta(velocity, theta[i], ee[i],\
        dbeta, dphi, uu, dphi_is_large[i])
        
      
  # the distribution of ejection angle is defined in coordinates (wpsi, wlambdaM) 
  # where wpsi is an angle between the jet main axis of symmetry and the direction of ejection
  # wlambdaM is a longitude in the plane perp# endicular to the jet's main axis
  # However, the factor 1/cos(psi) comes from the Jacobian of transformation 
  # (alphaM, betaM, u, psi, lambdaM) -> (alpha, beta, v, theta, lambda)
  # and here psi is an angle between the direction of ejection and the normal to surface

      #print(wpsi)
      #print(psi)
      #print(" wpsi = " + str(wpsi[i]), "psi = " + str(psi[i]))

      fac2 = distf.ejection_direction_distribution(var.source.ejection_angle_distr, wpsi[i],psi[i], lambdaM, \
        var.source.zeta, var.source.eta)

      fac2 = fac2 / np.cos(psi[i])

      #print("fac1 = " + str(fac1) + "fac2 = " + str(fac2) + "dd = " + str(ddphidtheta[i]))
                  
      tmpIntegrand = fac1 * fac2 / np.abs(ddphidtheta[i]) * rate[i]

      # to calculate the flux
      if(const.flux) :
        tmpIntegrand = tmpIntegrand * velocity * np.abs(np.cos(theta[i]))
      # endif

      #print("tmp = " +str(tmpIntegrand))
      
      Integrand = Integrand + tmpIntegrand

      if(tmpIntegrand  !=  tmpIntegrand) :
        f = open("PyPlumes/results/Integrand_number_density_output.txt","a")
        f.write("\n")
        f.write("\nNaN is obtained for an integrand value \n")
        f.write("factor related to ejection speed distribution: fac1 = " + str(fac1) + "\n")
        f.write("factor related to ejection direction distribution: fac2 = " + str(fac2) + "\n")
        f.write("the partial derivaive of delta phi by theta = " + str(ddphidtheta[i]) + "\n")
        f.write("theta = " + str(theta[i]) + "   psi = " + str(psi[i]) +"   lambdaM = "+ str(lambdaM) + "\n")
        f.write("dust production rate = " + str(rate[i]) +"\n") 
        f.close
      # endif
    # endif
    tmpIntegrand = 0.0
  # enddo
  #print(Integrand)
  return Integrand
# # end def Integrand_number_density  




# theta_geometry_hyperbola recieves vectors' absolute
# values and the angle between these vectors
# it's assumed that the point (0, 0, 0) is a focus of an hyperbola
# value of a semi major axis of the hyperbola, the points lay on is
# also an input parameter of the function
# the function returns theta that is the angle between radius-vector r
# and a tangent to the hyperbola in the point r
# the choice between two possible values of this angle is made
# in the way that movement from point rm to the point r along
# the hyperbola is possible
# (means: r and rm lay on the same quadrant in the CS in which
# the equation of the hyperbola is in canonical form)
# In general case there are two possible hyperbolae
# => two possible values of theta are to be investigated
# Returns array of 2 vallues of theta.
# If theta = large negative number,
# it means "no physically plausible solution can be found"
def theta_geometry_hyperbola(r0, rm0, vv, phi, a0,timeDependence):
  solved =  False 
  theta = np.array([-555.0, -555.0], dtype = float)
  dphi_is_large =  np.array([False,False]) 
  deltat = np.array([0.0,0.0], dtype =  float) 
  r2d = np.array([0.0,0.0])
  rm2d = np.array([0.0,0.0])
  shift = np.array([0.0,0.0])
                          
  r = r0 / rm0 ; rmoon = rm0 / rm0; 
  a = float(a0 / rm0)      				
  # define x-axis in the same direction as r-vector
  r2d[0] = r ; r2d[1] = 0.0    			
  
  # we don't have enough information to define the sign of r
  # vector in the right-handed coordinate system
  # but the value of theta that we are looking for are the same
  # in both cases
  
  rm2d[0] = rmoon * np.cos(phi) ; rm2d[1] = rmoon * np.sin(phi)
  
  # (x(1),y(1)) and (x(2),y(2)) are coordinates
  # of 2 possible position of the hyperbola's second focus
  x = np.array([0.0,0.0])
  y = np.array([0.0,0.0])
  x,y = ultilities.circle_intersection(rm2d[0], rm2d[1], 2.0 * a + rmoon, r2d[0], 2.0 * a + r)  	
  ee = np.array([0.0,0.0], dtype = float) 
    
  for i  in range(0, 2):
    # shift is coordinates of the ellipse's center
    # in the CS centered at the focus
    shift[0] = x[i] / 2.0 ; shift[1] = y[i]  / 2.0  	
    
    # coords of vector r in CS centered at center of the ellipse
    r2d = r2d - shift  
    # coords of vector rm in CS centered at center of the ellipse								
    rm2d = rm2d - shift      					
    
    # angle between major axis of the ellipse and the current x-axis
    angle = np.arctan(y[i] / x[i])
    # vector r in the CS with its center at the center of the ellipse
    # and the x-axis along major axis of the ellipse									
    r2d = ultilities.rot2d(r2d, -angle)  
    # vector rm in the CS with its center at the center of the ellipse
    # and the x-axis along major axis of the ellipse								
    rm2d = ultilities.rot2d(rm2d, -angle)  
    # vector shift in the CS with its center at the center of the ellipse
    # and the x-axis along major axis of the ellipse							
    shift = ultilities.rot2d(shift, -angle)


    # the hyperbola solves our problem only if r and rm lay
    # on the same branch and the trajectory doesn't intersect
    # the moon's surface (the particle doesn't pass the pericenter
    if (r2d[0] / rm2d[0] > 0.0 and r2d[1] / rm2d[1] > 0.0):  				
      c = float(0.50 * np.sqrt((x[i]**2) + (y[i]**2)))  
      #print(c**2)
      #print(a**2)
      b = float(np.sqrt(((c**2) - (a**2))))   						
      ee[i] = c / a
      one_plus_e = 1.0 + ee[i]
      one_minus_e = 1.0 - ee[i]
      one_minus_e2 = one_plus_e * one_minus_e
      aux = -a * one_minus_e2
      cosf1 = float((aux - 1.0) / ee[i])
      #print("aux is " + str(aux))
      #print("ee is " + str(ee))
      #print("cosf1 is " + str(cosf1))
      f1 = np.arccos(cosf1)
      
      f2 = f1 + phi
      
      theta[i] = const.halfpi - np.arctan((ee[i] * np.sin(f2)) / (1.0 + ee[i] * np.cos(f2)))
                  
      solved,discr = control(theta[i], ee[i], phi, vv, r0, rm0, dphi_is_large[i])

      if timeDependence == 1 :
        ean = 2.0 * np.arctanh(np.tan(f2/2.0) * np.sqrt(-one_minus_e / one_plus_e))
        eanm = 2.0 * np.arctanh(np.tan(f1/2.0) * np.sqrt(-one_minus_e / one_plus_e))
        tmp1 = a0 * np.sqrt(a0 / (const.gm)) * (ee(i) * np.sinh(eanm) - eanm)
        tmp2 = a0 * np.sqrt(a0 / const.gm) * (ee(i) * np.sinh(ean) - ean)
        deltat[i] = tmp2 - tmp1
      else:
        deltat[i] = 0.0
      # endif

      if (solved != 1 and  theta[i] > 0.0) :
        f = open("PyPlumes/results/theta_geometry_output.txt","a")
        f.write("\n")
        f.write("\nTheta was found with an insufficient accuracy of " + str(discr) + " from the geometry of hyperbola\n")
        f.write("for r = " + str(r0) + "\n")
        f.write("dphi = " + str(phi) + "\n")
        f.write("the obtained value of eccentricity = " + str(ee[i]) + "\n")
        f.write("and the obtained value of theta = " + str(theta[i]) + "\n")
        f.close() 
       
        if dphi_is_large[i] == 1 :
          f = open("PyPlumes/results/theta_geometry_output.txt","a")
          f.write("\nthe case of \Delta\phi > pi has been encoutered\n")
          f.close()
        # endif
      # endif
        
    else:
      theta[i] = -444.0
      ee[i] = -444.0
      dphi_is_large[i] =  0 
    # endif
    # we have changed the vectors r and rmoon we started from
    # so we need to go back to the beginning
    # to find the second value of theta
    r = r0 / rm0 ; rmoon = rm0 / rm0; a = a0 / rm0
    rm2d[0] = rmoon * np.cos(phi) ; rm2d[1] = rmoon * np.sin(phi)
    r2d[0] = r ; r2d[1] = 0.0
  
  return ee, theta, deltat, dphi_is_large
  # enddo

# # end def theta_geometry_hyperbola


# theta_geometry_ellipse recieves absolute values
# of two vectors and the angle between these vectors
# it's assumed that the point (0, 0, 0) is a focus of an ellipse
# value of a semi major axis of an ellipse, the points lay on
# is also an input parameter of the def
# the def returns theta that is the angle between
# radius-vector r and a tangent to the ellipse in the point r
# the choice between two possible values of this angle is made
# in the way that movement happens from point rm to the point r 
# In general case there are two possible ellipses
# => two possible values of theta are to be found
# Returns array of 2 vallues of theta.
# If theta = large negative number
# it means "no physically plausible solution can be found"

def theta_geometry_ellipse(r0, rm0, vv, phi, a0, timeDependence): 
  solved =  False 
  dphi_is_large =  np.array([False, False]) 
  theta = np.array([-888.0, -888.0])
  
  r = r0 / rm0  ; rmoon = rm0 / rm0
  a = a0 / rm0
  # define x-axis in the same direction as rmoonvector
  r2d = np.array([0.0,0.0])
  r2d[0] = r; r2d[1] = 0.0    			
  # we don't have enough information to define the sign of r
  # vector in the right-handed coordinate system
  # but the value of theta that we are looking for are the same
  # in both cases
  rm2d = np.array([0.0,0.0])
  rm2d[0] = rmoon * np.cos(phi) ; rm2d[1] = rmoon * np.sin(phi)
  
  # (x(1),y(1)) and (x(2),y(2)) are coordinates of 2 possible
  # position of the ellsipse's second focus
  x, y = ultilities.circle_intersection(rm2d[0], rm2d[1], 2.0 * a - rmoon, r2d[0], 2.0 * a - r) 
  # distance between the foci of the ellipse
  cc = np.sqrt((x**2) + (y**2))
  ee = np.array([0.0,0.0])
  deltat = np.array([0.0,0.0])  
  
  for i  in range(0, 2):
    if cc[i] == cc[i]  :
      ee[i] = cc[i] / 2.0 / a

      one_plus_e = 1.0 + ee[i]
      one_minus_e = 1.0 - ee[i]
      one_minus_e2 = one_plus_e * one_minus_e
      cosf1 = (-x[i] * rm2d[0] - y[i] * rm2d[1]) / rmoon / cc[i]
      f1 = np.arccos(cosf1)
      f2 = f1 + phi
      
      if(r < rmoon) :
        rtest = a * one_minus_e2 / (1.0 + ee[i] * np.cos(f2))
        # if r and rm are both located after apocenter
        if( np.abs(1.0 - rtest / r ) > 1*(10**-6)) :  			
          f2 = (const.twopi - f1) + phi
          rtest = a * one_minus_e2 / (1.0 + ee[i] * np.cos(f2))
          # rm is before apocenter, the case of large dphi encountered	
          # phi is the angle between vectors r and rm.
          # it can be that between r and rm the particle
          # traveled the arc of 2pi - phi
          # if the direction of movement along the ellipse
          # is chosen incorrect rtest  !=  r
          if( np.abs(1.0 - rtest / r ) >  1*(10**-6)) :  		
            f2 = f1 + (2.0 * const.pi - phi)
            dphi_is_large[i] =  1 
            rtest = a * one_minus_e2 / (1.0 + ee[i] * np.cos(f2))
            # rm is after apocenter, the case of large dphi encountered
            if( np.abs(1.0 - rtest / r ) >  1*(10**-6)) :  	
              f1 = const.twopi -f1
              f2 = f1 + (2.0 * const.pi - phi)
              dphi_is_large[i] =  1
            # endif
          else:
            f1 = const.twopi - f1
          # endif
        # endif
        if(f2 > const.twopi):
           f2 = f2 - const.twopi
      else:  					
        rtest = a * one_minus_e2 / (1.0 + ee[i] * np.cos(f2))
        if( np.abs(1.0 - rtest / r ) >  1*(10**-6)) :
          f2 = f1 + (2.0 * const.pi - phi)
          dphi_is_large[i] =  1
        # endif
      # endif

      # No ejection downward even if rM > r
      if(f1 < const.pi - 1*(10**-4)) :  		
        theta[i] = const.halfpi - np.arctan((ee[i] * np.sin(f2)) \
                    / (1.0 + ee[i]* np.cos(f2)))
        
        solved, discr = control(theta[i], ee[i], phi, vv, r0, rm0, \
              dphi_is_large[i])
      
        if( timeDependence == True) :
          ean = 2.0 * np.arctan(np.tan(f2/2.0) \
              * np.sqrt(one_minus_e / one_plus_e))
          if(ean < 0.0):
            ean = const.twopi + ean

          eanm = 2.0 * np.arctan(np.tan(f1/2.0) \
              * np.sqrt(one_minus_e / one_plus_e))
          tmp1 = a0 * np.sqrt(a0 / const.gm) * (eanm - ee[i] * np.sin(eanm))
          tmp2 = a0 * np.sqrt(a0 / const.gm) * (ean - ee[i]* np.sin(ean))
          deltat[i] = tmp2 - tmp1
        else:
          deltat[i] = 0.0
        # endif
      else:
        theta = np.array([-777.0,-777.0])
      # endif
                  
      if( solved != 1  and  theta[i] > 0.0) :
        f = open("PyPlumes/results/theta_geometry_output.txt","a")
        f.write("\n")
        f.write("\nTheta was found with an insufficient accuracy of " + str(discr) + " from the geometry of ellipse\n")
        f.write("for r = " + str(r0) + "\n")
        f.write("rt = " + str(rtest * rm0) + "\n")
        f.write("f1 = " + str(f1) + "\n") 
        f.write("f2 = " + str(f2) + "\n")
        f.write("dphi = " + str(phi) + "\n")
        f.write("the obtained value of eccentricity = " + str(ee[i]) + "\n") 
        f.write("and the obtained value of theta = " + str(theta[i]) + "\n")
        f.close() 
       
        if dphi_is_large[i] == 1 :
          f = open("PyPlumes/results/theta_geometry_output.txt","a")
          f.write("\nthe case of \Delta\phi > pi has been encoutered\n")
          f.close()
    
        # endif
      # endif
    else:
      theta[i] = -888.0
      ee[i] = -888.0
      dphi_is_large[i] =  False 
    # endif
  return ee, theta, deltat, dphi_is_large
  # enddo


# # end def theta_geometry_ellipse





# def control tests if the obtained value of theta is correct
# the criterion is: using  the obtained value of theta
# one gets the same value of dphi which was used to calculate the theta
# also the eccentricity values are compared
# input parameters: r0 and rm0 - lengths of two vectors - positions on the orbit,
# vv - speed at position r0, phi - angle between r0 and rm0 used
# in def theta_geometry_...to calculate theta
# theta is the angle between r0 and velosity at the position r0
# ee is eccentricity obtained in theta_geometry...
# dphi_is_large tells if it is the case when the particle traveled
# from rm to r over an arc of 2pi - dphi
# the def control uses formulae for energy
# and angular momentum to obtain values of ee1 and dphi
# they must differ in less : eps, in this case the def returns TRUE
# otherwise it returns FALSE in the variable solved
# it also returns the estimated accuracy stored in the variable discr

def control(theta, ee, phi, vv, r0, rm0, dphi_is_large):
 
  Ekep = vv * vv / 2.0 - const.gm / r0
  hh = r0 * vv * np.sin(theta)
  hh2 = hh * hh
  solved = False
  eps = 1 * (10**-4)

#  eccentricity (eq 31)
  ee1 = float(np.sqrt(1.0 + 2.0 * Ekep * (hh / const.gm) * (hh / const.gm)))
  if(np.abs(ee1 - ee) > eps) :
    f = open("PyPlumes/results/theta_geometry_output.txt","a")
    f.write("\neccentricity is incorrect: " + str(ee) + " instead of " + str(ee1) + "\n")
    f.close
    # endif

#   delta phi (equation 32)
  cosp = (hh2 / r0 / const.gm - 1.0) / ee1
  cospm = (hh2 / rm0 / const.gm - 1.0) / ee1
  #######################################
  if(cosp > 1.0) :
    cosp = 1.0
    f = open("PyPlumes/results/theta_geometry_output.txt","a")
    f.write("\ncos(phi) > 1 obtained, corrections applied\n")
    f.close 
  
  # endif
  if(cosp < -1.0) :
    cosp = -1.0
    f = open("PyPlumes/results/theta_geometry_output.txt","a")
    f.write("\ncos(phi) < -1 obtained, corrections applied\n")
    f.close
  # endif
  
  if(cospm > 1.0) :
    cospm = 1.0
    f = open("PyPlumes/results/theta_geometry_output.txt","a")
    f.write("\ncos(phiM) > 1 obtained, corrections applied\n")
    f.close 
  # endif

  if(cospm < -1.0) :
     cospm = -1.0
     f = open("PyPlumes/results/theta_geometry_output.txt","a")
     f.write("\ncos(phiM) < -1 obtained, corrections applied\n")
     f.close 
    
  # endif
  ########################################
  phi1 = np.arccos(cosp)
  phi1m = np.arccos(cospm)
  if(theta < const.halfpi) :
    dphi = phi1 - phi1m
  else:
    dphi = (2.0 * const.pi - phi1) - phi1m
  # endif
  
  if(dphi_is_large) :
    discr = np.abs(2.0 * const.pi - dphi - phi) / np.abs(phi)
    if discr < eps :
      solved == True
  else:
    discr = np.abs(dphi - phi) / np.abs(phi)
    if discr < eps: 
      solved == True
  
  return solved, discr
  # endif
  
# # end def control
  
  



# end module twobody_fun
