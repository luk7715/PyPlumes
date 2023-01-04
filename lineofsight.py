import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys
import inspect
from scipy import integrate as integrate

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import variables as var
import integrator
import ultilities as ult
import input
import output
import distributions as distf
import gu
import const

Vrad = 4.5* 10** (-5)
Hrad = 4.5* 10** (-5)
Vpix = 64
Hpix = 64
ni = 20
sampdist_large = 3**2 #meters


#returns an array of spacecraft positions placed along the lines of sight
#of fictive pixels of a fictive CCD frame
def line_of_sight(source): 
  
    Hcamscale = Hrad / float(Hpix)
    Vcamscale = Vrad / float(Vpix) # rad / pixel in horisontal and vertical direction

    scpos = np.array([0.0,0.0,0.0])
    scpos[0] = 6.282*(10**10) * const.sqrt2 * 3.0
    scpos[1] = 6.282*(10**10) * const.sqrt2 * 3.0
    scpos[2] = 6.282*(10**10) * 1.120
    scpos = scpos
    
 
    xcam = np.array([0.0,0.0,1.0])
    ycam = np.array([ -1.0/const.sqrt2 , 1.0/const.sqrt2 ,0.0])
    zcam = np.array([ -1.0/const.sqrt2 ,-1.0/const.sqrt2 , 0.0])
    # distance from spacecraft to the "center of the plume", meters
    r0 = ult.norma3d(scpos-var.source.rrM)
    
    # moon angular radius as appears from the point of observation
    arcrm = const.rm / ult.norma3d(scpos)

    # direction towards moon center in camera CS
    mooncenterdir = np.array([(-scpos[2]/const.rm), (0.0), (scpos[0] * const.sqrt2 / const.rm) ])
    mooncenterdir = mooncenterdir / ult.norma3d(mooncenterdir)
    pixdircam = np.array([0.0,0.0,0.0])
    pixdir = np.array([0.0,0.0,0.0])
    
    for i in range(-Hpix, Hpix):
            # pixdircam is a unit vector in camera's CS pointing in the direction which is pictured in the pixel with coords (i,ii)
            # (0,0,1) would be a direction of the center
            # of the matrix where 4 central pixels meet
            
            if(i < 0):
                 pixdircam[0] = (float(i)+0.50) * Hcamscale
            if(i > 0):
                 pixdircam[0] = (float(i)-0.50) * Hcamscale	
                        
            for ii in range (-Vpix, Vpix):

                if(ii < 0):
                     pixdircam[1]= (float(ii)+0.50) * Vcamscale
                if(ii > 0):
                     pixdircam[1] = (float(ii)-0.50) * Vcamscale
                
                pixdircam[2] = np.sqrt(1.0 - pixdircam[0]**2 - pixdircam[1]**2)
                # angle between 2 unit vectors
                arctmp_pre = pixdircam * mooncenterdir
                #print(arctmp_pre)
               
                arctmp = np.arccos(np.sum(arctmp_pre))
                #print(str("arctmp is " )+ str(arctmp))
                #print(str("arcrm is " )+ str(arcrm))
                # exclude disc of the moon and center cross of the array elements
                # of which don't correspond to fictive or real pixels
                if(arctmp >= arcrm and i != 0 and ii != 0):
                    #print("calculate")
                    
                    # pixdir is pixdircam in the moon-centered CS
                    pixdir[0] = np.sum(pixdircam * np.array([xcam[0], ycam[0], zcam[0]]))
                    pixdir[1] = np.sum(pixdircam * np.array([xcam[1], ycam[1], zcam[1]]))
                    pixdir[2] = np.sum(pixdircam * np.array([xcam[2], ycam[2], zcam[1]]))

                    distance, K1, K2= ult.dist_between_2lines(scpos, pixdir, var.source.rrM, var.source.symmetry_axis)

                    for iii in range(-ni,ni):
                        rvector = K1 + pixdir * iii * sampdist_large
                        r = ult.norma3d(rvector)
                        r_scaled = r/const.rm
                        alpha = np.arccos(rvector[2]/r)
                        beta = ult.myatan1(rvector[0],rvector[1])
                        #with the given maximum ejection velocity the radial distance
                        #which the dust particle can achieve is restricted
                        #also the density inside the moon should not be computed
                        compute1 = bool(0)
                        #print(r)
                        if( r < 6.282*(10**13) and r > const.rm):
                            compute1 = True 
                            #print("cal")

                        var.pointS[i,ii,iii]= var.point(r, r_scaled, alpha, beta, rvector, compute1)
                        #print(var.pointS[i,ii,iii].compute)
                 
                
                else:
                # if the line of sight points to the moon's disk or the pixels with a 0-index
                    for iii in range(-ni,ni):
                        compute = bool(0)
                        compute = False
                        var.pointS[i,ii,iii]= var.point(0, 0, 0, 0, 0, compute)
                    

    return var.pointS