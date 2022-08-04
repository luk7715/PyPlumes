import numpy as np
import pandas as pd
import const
import variables as var 
import gu
import help 

# input the parameters of the sources from the file with the name
# stored in the variable fname
def read_sources_params(fname, Ns):

  data = pd.read_csv(fname, nrows = Ns)
  data_in_arrays = data.values
  
  var.source.alphaM = data_in_arrays[:,0] ; var.source.betaM = data_in_arrays[:,1]
  var.source.zeta = data_in_arrays[:,2] ; var.source.eta = data_in_arrays[:,3]
  var.source.production_fun = data_in_arrays[:,4] ; var.source.production_rate = data_in_arrays[:,5]
  var.source.ud_shape = data_in_arrays[:,6] ; var.source.ud_umin = data_in_arrays[:,7]
  var.source.ud_umax = data_in_arrays[:,8] ; var.source.ejection_angle_distr = data_in_arrays[:,9]
  var.source.sd = data_in_arrays[:,10]

  var.source.r = const.rm
  var.source.alphaM = var.source.alphaM * const.deg2rad
  var.source.betaM = var.source.betaM * const.deg2rad
  var.source.alphaM = const.halfpi - var.source.alphaM

  var.source.zeta = var.source.zeta *const.deg2rad
  var.source.eta = var.source.eta * const.deg2rad

  var.source.rrM[0] = const.rm *np.sin(var.source.alphaM) * np.cos(var.source.betaM)
  var.source.rrM[1] = const.rm *np.sin(var.source.alphaM) * np.sin(var.source.betaM)
  var.source.rrM[2] = const.rm *np.cos(var.source.alphaM)

  var.source.is_jet = True

  ui, Si = gu.Gu_integral(var.source.sd)
  var.source.ui = ui
  var.source.Gu_precalc = Si 

  axis = jet_direction(var.source.betaM, var.source.zeta, var.source.eta,var.source.rrM)
  var.source.symmetry_axis = axis 
          				
# end function read_sources_params



# obtains Cartesian coordinates of a unit vector 
# in the moon-centered coordinate system
# aligned with the ejection symmetry axis
# the source's position (rrM - Cartesian coordinates, betaM - eastern longitude),
# the jet's zenith angle (zeta) and azimuth (eta) are known
def jet_direction(betaM, zeta, eta, rrM):

  rtmp = rrM / const.rm
  tmpang = 3.0 * const.halfpi-betaM
  if(zeta != 0.0): 
    xout, yout, zout = help.eulrot(0.0, 0.0, tmpang, rtmp[0], rtmp[1], rtmp[2], 0)
    rtmp[0] = xout ; rtmp[1] = yout ; rtmp[2] = zout
    
    xj[0] = 0.0
    xj[1] = 1.0 * np.sign(rtmp[2]) * np.abs(rtmp[2])
    xj[2] = -1.0 * np.sign(rtmp[1]) * np.abs(rtmp[1])
    xj = xj/help.norma3d(xj)

    yj = help.vector_product(rtmp,xj)
        
    jetdir = np.sin(zeta) * np.cos(eta) * xj - np.sin(zeta) * np.sin(eta) * yj \
        + np.cos(zeta) * rtmp

    jetdir = jetdir / help.norma3d(jetdir)
    
    xout, yout, zout = help.eulrot(0.0, 0.0, tmpang, rtmp[0], rtmp[1], rtmp[2],1)

    rtmp[0]= xout ; rtmp[1] = yout ; rtmp[2] = zout
    
    xout, yout, zout = help.eulrot(0.0, 0.0, tmpang, jetdir[0], jetdir[1], jetdir[2], 1)

    jetdir[0] = xout ; jetdir[1] = yout ; jetdir[2] = zout

  else: 
    jetdir = rtmp
  
  return jetdir