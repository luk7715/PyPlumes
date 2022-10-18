import const
import variables as var
import numpy as np


def result_out(density,nt):
  '''
  This function is used to write output results to a text file.
  
  Args:
   -density: The output array of density values
   -nt: Number of points 
  '''

  f = open("PyPlumes/results/twobody_model_result.txt","a") 
  for i  in range(0, nt):
    f.write(str(density[i,:]) + " " + str(density[i,:] )+ " " + str(var.point.r) + \
       " " + str((const.halfpi - var.point.alpha[i]) * const.rad2deg) \
        + " " + str(var.point.beta[i]*const.rad2deg))
    f.write("\n")
  f.close()





def cassini_flyby_out(density, ttab, bg, nt):
  '''
  This is specifically written for writing the output of the enceladus example. See file "./testing/enceladus_test.py"
  '''
  
  f = open("PyPlumes/results/E2_profile.dat","w")
  for i  in range(0, nt):
    densitysum = np.sum (density[i,:]) + bg  
    f.write(str(ttab[i]) + " " + str(densitysum) + "\n")
  f.close 
# # end def cassini_flyby_out  


def surface_deposition_out(num, deposition, nt, dphi):
  '''
  This is specifically written for writing the output of the europa example. See file "./testing/europa_test.py"
  '''
  
  match num:
    case 1:
      fname = 'PyPlumes/results/narrow_jet_shallow_sd.dat'
    case 2:
      fname = 'PyPlumes/results/diffuse_source_steep_sd.dat'
    case 3:
      fname = 'PyPlumes/results/diffuse_source_shallow_sd.dat'
    case 4:
      fname = 'PyPlumes/results/narrow_jet_steep_sd.dat'
   
  f = open(fname,"w") 
  for i  in range(0, nt):
    f.write(str(dphi[i] *const.rm * 10**(-3)) + " " + str(deposition[i])+ " " )
    f.write("\n")
  f.close()
			



