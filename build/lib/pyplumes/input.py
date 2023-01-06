import numpy as np
import pandas as pd
import const
import variables as var 
import gu
import ultilities as ultilities 



def read_spacecraft_coordinates(fname,nt):
  '''
  This function is the template for inputing spacecraft coordinates from files.
  In the file fname the coordinates should be written as:
  radial distance from the moon center [m], latitude [deg], eastern longitude [deg].
  The class var.point will be constructed using the data input.

  Args:
   -fname: Path to the desired file, string
   -nt: Number of points to be calculated, integer
  '''
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



def read_Cassini_E2(nt):
  '''
  This function is specifically written for reading in the Cassini E2 flyby data based 
  on the template above. The file was prepared in advance using SPICE software
  
  Args:
   -nt: Number of points to be calculated, integer
  '''
  data = np.loadtxt("PyPlumes/input_data_files/Cassini_E2_flyby.dat", max_rows= nt)
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




def read_sources_params(fname, Ns):
  '''
  This function is the template for inputing the parameters of the sources from files.
  The class var.source will be constructed using the data input.

  Args:
   -fname: Path to the desired file, string
   -Ns: Number of sources, integer
  '''

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



def get_europa_input(Ns, nt):
  '''
  This function is specifically written for getting the input data for the Europa
  surface depositions calculation based on the templates above.

  Args:
   -nt: Numer of points to be calculated, integer
   -Ns: Number of sources, integer
  '''
  ead_choice = np.array([1, 3, 3, 1])
  sd_choice = np.array([2, 3, 2, 3])
  
  # define 4 sources with the same coordinates and verticle axis
  # of symmetry but different size- and ejection direction distributions
  var.source.alphaM= const.halfpi; var.source.betaM = 0.0; var.source.zeta=0.0
  var.source.eta= 0.0; var.source.production_fun= 0
  var.source.production_rate = 1*(10**14)
  var.source.ud_shape = 1
  var.source.ud_umin = 0.0
  var.source.ud_umax = 500.0
  var.ud.ud_shape = 1 ; var.ud.umin = 0.0
  var.ud.umax = 500.0
  var.source.r = const.rm
  var.source.ejection_angle_distr = ead_choice
  var.source.sd = sd_choice
  var.source.rrM = np.array([0.0,0.0,0.0])
  var.source.rrM[0] = const.rm *np.sin(var.source.alphaM) * np.cos(var.source.betaM)
  var.source.rrM[1] = const.rm *np.sin(var.source.alphaM) * np.sin(var.source.betaM)
  var.source.rrM[2] = const.rm *np.cos(var.source.alphaM)
  var.source.is_jet = True

  
  ui = np.zeros([const.GRN, Ns])
  Si = np.zeros([const.GRN, Ns])
  for i  in range(0, Ns):
    ui[:,i], Si[:,i] = gu.Gu_integral(var.source.sd[i])
    #print(ui[:,i])

  #print("ui is " + str(ui))
  #print(Si)

  var.source.ui = ui
  var.source.Gu_precalc = Si 

  axis = jet_direction(var.source.betaM, var.source.zeta, var.source.eta,var.source.rrM)
  var.source.symmetry_axis = axis 
  
  # the points are equidistantly placed on a 10deg arc
  # having the source at one of the arc's # ends
  dphi = np.zeros(nt)

  for i in range (0,nt):
    dphi[i] = 2.0 * float(i+1) / float(nt) *const.deg2rad
 
  var.point.r = const.rm
  var.point.alpha = const.halfpi 
  var.point.beta = dphi
  var.point.rvector = np.zeros([nt,3])   

  #print (var.point.r * np.sin(var.point.alpha) * np.cos(var.point.beta))
  for i  in range(0, nt):
    var.point.rvector[i,0] = var.point.r * np.sin(var.point.alpha) * np.cos(var.point.beta[i])
    var.point.rvector[i,1] = var.point.r * np.sin(var.point.alpha) * np.sin(var.point.beta[i])
    var.point.rvector[i,2] = var.point.r * np.cos(var.point.alpha) 
  
  var.point.r_scaled = 1.0
  var.point.compute = True

  

  return dphi, var.point, var.source

# end subroutine get_europa_input

def europascatter_input2(nt):
  '''
  This function is specifically written for getting the input data for the Europa
  surface depositions calculation based on the templates above.

  Args:
   -nt: Numer of points to be calculated, integer
   -Ns: Number of sources, integer
  '''
 
  ky_production = np.array([1.63043478e+21, 1.30434783e+19, 1.63043478e+18,1.30434783e+16, 1.63043478e+15])
  # define 4 sources with the same coordinates and verticle axis
  # of symmetry but different size- and ejection direction distributions
  var.source.alphaM= const.halfpi; var.source.betaM = 0.0; var.source.zeta=0.0
  var.source.eta= 0.0; var.source.production_fun= 0
  var.source.production_rate = ky_production[4]
  var.source.ud_shape = 1
  var.source.ud_umin = 0.0
  var.source.ud_umax = 1500.0
  var.ud.ud_shape = 1 ; var.ud.umin = 0.0
  var.ud.umax = 1500.0
  var.source.r = const.rm
  var.source.ejection_angle_distr = 1
  var.source.sd = 4
  var.source.rrM = np.array([0.0,0.0,0.0])
  var.source.rrM[0] = const.rm *np.sin(var.source.alphaM) * np.cos(var.source.betaM)
  var.source.rrM[1] = const.rm *np.sin(var.source.alphaM) * np.sin(var.source.betaM)
  var.source.rrM[2] = const.rm *np.cos(var.source.alphaM)
  var.source.is_jet = True
  var.source.symmetry_axis = np.array([0.0,0.0,0.0])

  

  ui, Si = gu.Gu_integral(var.source.sd)
 
  #print("ui is " + str(ui))
  #print(Si)

  var.source.ui = ui
  var.source.Gu_precalc = Si 
  
  axis = jet_direction(var.source.betaM, var.source.zeta, var.source.eta,var.source.rrM)
  var.source.symmetry_axis = axis 
  
  # the points are equidistantly placed on a 10deg arc
  # having the source at one of the arc's # ends
  dphi = np.zeros(nt)

  for i in range (0,nt):
    dphi[i] = 2.0 * float(i+1) / float(nt) *const.deg2rad

  

  return dphi, var.source

# end subroutine get_europa_input

def jet_direction(betaM, zeta, eta, rrM):
  '''
  This function is used to obtain the Cartesian coordinats of a unit vector
  in the moon-centered coordinate system aligned with the ejection symmetry axis.
  The source position data is used for this calculation.

  Args:
   -betaM: Eastern Longitude of the jet
   -zeta: The zenith angle of the jet
   -eta: The azimuth of the jet
   -rrM: Cartesian coordinates of the jet
  Return:
   -jetdir:  (3d-vector) unit vector in moon-centered coordinate system pointing to the direction of 
            the axis around which ejection is symmetrical
  Return Type:
   -Array
  '''

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
