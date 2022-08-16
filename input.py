import numpy as np
import pandas as pd
import const
import variables as var 
import gu
import ultilities as ultilities 



# A template for inputing the spacecraft coordinates from the file with the name 
# stored in the variable fname
# in the file fname the coordinates are written as follows
# radial distance from the moon center [m], latitude [deg], eastern longitude [deg]
def read_spacecraft_coordinates(fname,nt):
  data = np.loadtxt(fname, max_rows= nt)
  data_in_arrays = data
  #print(data_in_arrays)


  var.point.r = data_in_arrays[:,0] ; var.point.alpha = data_in_arrays[:,1]
  var.point.beta = data_in_arrays[:,2]

  var.point.alpha = var.point.alpha * const.deg2rad
  var.point.alpha = const.halfpi - var.point.alpha 
  var.point.beta = var.point.beta *const.deg2rad

  var.point.rvector = np.zeros([nt,3])
  for i in range (0,nt):
    var.point.rvector[i,0] = var.point.r[i] * np.sin(var.point.alpha[i]) * np.cos(var.point.beta[i])
    var.point.rvector[i,1] = var.point.r[i] * np.sin(var.point.alpha[i]) * np.sin(var.point.beta[i])
    var.point.rvector[i,2] = var.point.r[i] * np.cos(var.point.alpha[i]) 

  var.point.r_scaled = var.point.r/const.rm
  var.point.compute = True 

  return var.point
# end subroutine read_spacecraft_coordinates


##Specifically for read in the Cassini E2 flyby data file based on the template above
#the file was prepared in advance using SPICE software
def read_Cassini_E2(nt):
  data = np.loadtxt("PyPlumes/input_data_files/Cassini_flyby_test.dat", max_rows= nt)
  data_in_arrays = data
  #print(data_in_arrays)

  ttab = np.zeros([nt,1])
  ttab = data_in_arrays[:,0]

  var.point.r = data_in_arrays[:,1] ; var.point.alpha = data_in_arrays[:,2]
  var.point.beta = data_in_arrays[:,3]

  var.point.alpha = var.point.alpha * const.deg2rad
  var.point.alpha = const.halfpi - var.point.alpha 
  var.point.beta = var.point.beta *const.deg2rad
  #print("test")
  #print(var.point.r[0] * np.sin(var.point.alpha[0]) * np.cos(var.point.beta[0]))

  var.point.rvector = np.zeros([nt,3])
  for i in range(0, nt):
   var.point.rvector[i,0] = var.point.r[i] * np.sin(var.point.alpha[i]) * np.cos(var.point.beta[i])
   var.point.rvector[i,1] = var.point.r[i] * np.sin(var.point.alpha[i]) * np.sin(var.point.beta[i])
   var.point.rvector[i,2] = var.point.r[i] * np.cos(var.point.alpha[i]) 

  var.point.r_scaled = var.point.r/const.rm
  var.point.compute = True 

  return ttab, var.point


# input the parameters of the sources from the file with the name
# stored in the variable fname
def read_sources_params(fname, Ns):

  data = np.loadtxt(fname, max_rows= Ns)
  data_in_arrays = data
  #print(data_in_arrays)

  var.source.alphaM = data_in_arrays[0] ; var.source.betaM = data_in_arrays[1]
  var.source.zeta = data_in_arrays[2] ; var.source.eta = data_in_arrays[3]
  var.source.production_fun = data_in_arrays[4] ; var.source.production_rate = data_in_arrays[5]
  var.source.ud_shape = data_in_arrays[6] ; var.source.ud_umin = data_in_arrays[7]
  var.source.ud_umax = data_in_arrays[8] ; var.source.ejection_angle_distr = data_in_arrays[9]
  var.source.sd = data_in_arrays[10]

  var.ud.ud_shape = data_in_arrays[6] ; var.ud.umin = data_in_arrays[7]
  var.ud.umax = data_in_arrays[8] 
  #print("check input")
  #print (var.source.alphaM)

  var.source.r = const.rm
  var.source.alphaM = var.source.alphaM * const.deg2rad
  #print("in the function")
  #print (var.source.alphaM)
  var.source.betaM = var.source.betaM * const.deg2rad
  var.source.alphaM = const.halfpi - var.source.alphaM
  #print("check output")
  #print (var.source.alphaM)

  var.source.zeta = var.source.zeta *const.deg2rad
  var.source.eta = var.source.eta * const.deg2rad

  var.source.rrM = np.array([0.0,0.0,0.0])
  var.source.rrM[0] = const.rm *np.sin(var.source.alphaM) * np.cos(var.source.betaM)
  var.source.rrM[1] = const.rm *np.sin(var.source.alphaM) * np.sin(var.source.betaM)
  var.source.rrM[2] = const.rm *np.cos(var.source.alphaM)

  var.source.is_jet = True

  ui, Si = gu.Gu_integral(var.source.sd)
  var.source.ui = ui
  var.source.Gu_precalc = Si 

  axis = jet_direction(var.source.betaM, var.source.zeta, var.source.eta,var.source.rrM)
  var.source.symmetry_axis = axis 
  return var.source
          				
# end function read_sources_params



# obtains Cartesian coordinates of a unit vector 
# in the moon-centered coordinate system
# aligned with the ejection symmetry axis
# the source's position (rrM - Cartesian coordinates, betaM - eastern longitude),
# the jet's zenith angle (zeta) and azimuth (eta) are known
def jet_direction(betaM, zeta, eta, rrM):

  rtmp = rrM / const.rm
  tmpang = 3.0 * const.halfpi-betaM
  xj = np.array([0.0,0.0,0.0])
  jetdir = np.array([0.0,0.0,0.0])

  if(zeta != 0.0): 
    xout, yout, zout = ultilities.eulrot(0.0, 0.0, tmpang, rtmp[0], rtmp[1], rtmp[2], 0)
    rtmp[0] = xout ; rtmp[1] = yout ; rtmp[2] = zout
    
    xj[0] = 0.0
    xj[1] = 1.0 * np.sign(rtmp[2]) * np.abs(rtmp[2])
    xj[2] = -1.0 * np.sign(rtmp[1]) * np.abs(rtmp[1])
    xj = xj/ultilities.norma3d(xj)

    yj = ultilities.vector_product(rtmp,xj)
        
    jetdir = np.sin(zeta) * np.cos(eta) * xj - np.sin(zeta) * np.sin(eta) * yj \
        + np.cos(zeta) * rtmp

    jetdir = jetdir / ultilities.norma3d(jetdir)
    
    xout, yout, zout = ultilities.eulrot(0.0, 0.0, tmpang, rtmp[0], rtmp[1], rtmp[2],1)

    rtmp[0]= xout ; rtmp[1] = yout ; rtmp[2] = zout
    
    xout, yout, zout = ultilities.eulrot(0.0, 0.0, tmpang, jetdir[0], jetdir[1], jetdir[2], 1)

    jetdir[0] = xout ; jetdir[1] = yout ; jetdir[2] = zout

  else: 
    jetdir = rtmp
  
  return jetdir
