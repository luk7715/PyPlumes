# This file is a part of DUDI, the Fortran-95 implementation 
# of the two-body model for dust dynamics
# Version 1.0.0
# This is free software. You can use and redistribute it 
# under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
# If you do, please cite the following paper
# Anastasiia Ershova and Jürgen Schmidt, 
# Two-body model for the spatial distribution of dust ejected from
# an atmosphereless body, 2021, A&A, 650, A186 
# File: image_construction.f95
# Description: Additional parameters and subroutines used to construct an image

# Author: Anastasiia Ershova
# E-mail: vveyzaa@gmail.com

import const
import variables as var
import ultilities as ult
import numpy as np

## Define the variables 
Vrad = 4.5*(10**-2)
Hrad = 4.5*(10**-2)
Vpix = 64
Hpix = 64
ni = 20
sampdist_large = 3*10**2 	# meters

# returns an array of spacecraft positions placed alont the lines of sight
# of fictive pixels of a fictive CCD frame
def line_of_sight(points, source):
  #integer i, ii, iii, nn
  #real(8) euler_rot(3), xout, yout, zout
  #type(position_in_space), intent(out) :: points(-Hpix:Hpix,-Vpix:Vpix,-ni:ni)
  #type(source_properties), intent(in) :: source
  #real(8) Hcamscale, Vcamscale  	# rad / pixel in horisontal and vertical direction
  #real(8) trmat(3,3), pixdir(3), pixdircam(3), distance, K1(3), K2(3), scpos(3)
  #real(8) r0, arctmp, arcrm, mooncenterdir(3), tmp, xcam(3), ycam(3), zcam(3)
  
  Hcamscale = Hrad / float(Hpix)
  Vcamscale = Vrad / float(Vpix)
  scpos[0] = const.rm * np.sqrt(2.0) * 3.0
  scpos[1] = const.rm * np.sqrt(2.0) * 3.0
  scpos[2] = const.rm * 1.120
  scpos = scpos
  
  xcam = [00, 00, 10]
  ycam = [-10 / np.sqrt(2.0), 10 / np.sqrt(2.0), 00]
  zcam = [-10 / np.sqrt(2.0), -10 / np.sqrt(2.0), 00]
  # distance from Cassini to the "center of the plume", meters
  r0 = norma3d(scpos-source%rrM)
  
  # moon angular radius as appears from the point of observation
  arcrm = rm / norma3d(scpos)
  # direction towards moon center in camera CS
  mooncenterdir = [-scpos(3)/rm, 00, scpos(1) * sqrt20 / rm ]
  mooncenterdir = mooncenterdir / norma3d(mooncenterdir)
  
  
  for i  in range(-Hpix, Hpix):
      # pixdircam is a unit vector in camera's CS pointing
      # in the direction which is pictured in the pixel with coords (i,ii)
      # (0,0,1) would be a direction of the center
      # of the matrix where 4 central pixels meet
      
      if(i < 0) pixdircam(1) = (dble(i)+0.50) * Hcamscale
      if(i > 0) pixdircam(1) = (dble(i)-0.50) * Hcamscale  
            
    for ii  in range(-Vpix, Vpix):

      if(ii < 0) pixdircam(2) = (dble(ii)+0.50) * Vcamscale
      if(ii > 0) pixdircam(2) = (dble(ii)-0.50) * Vcamscale
      
      pixdircam(3) = sqrt(10 - pixdircam(1)**2 - pixdircam(2)**2)
      # angle between 2 unit vectors
      arctmp = acos(sum(pixdircam * mooncenterdir))
      # exclude disc of the moon and center cross of the array elements
      # of which don't correspond to fictive or #real pixels
      if(arctmp >= arcrm  and  i  !=  0  and  ii  !=  0) :
        
        # pixdir is pixdircam in the moon-centered CS
        pixdir(1) = sum(pixdircam * [xcam(1), ycam(1), zcam(1)])
        pixdir(2) = sum(pixdircam * [xcam(2), ycam(2), zcam(2)])
        pixdir(3) = sum(pixdircam * [xcam(3), ycam(3), zcam(3)])
        call dist_between_2lines(distance, K1, K2, scpos, pixdir, \
                                  source%rrM, source%symmetry_axis)
        
        forall(iii = -ni:ni)
          points(i,ii,iii)%rvector = K1 + pixdir * iii * sampdist_large
        # endforall
                  # if the line of sight is special
        for iii  in range(-ni, ni):
          points(i,ii,iii)%r = norma3d(points(i,ii,iii)%rvector)
          points(i,ii,iii)%r_scaled = points(i,ii,iii)%r / rm
          points(i,ii,iii)%alpha = acos(points(i,ii,iii)%rvector(3) \
                                        / points(i,ii,iii)%r)
          points(i,ii,iii)%beta = myatan1(points(i,ii,iii)%rvector(1), \
                                    points(i,ii,iii)%rvector(2))
          # with the given maximum ejection velocity the radial distance
          # which the dust particle can achieve is restricted
          # also the density inside the moon should not be computed
          points(i,ii,iii)%compute = points(i,ii,iii)%r < 1.98d6 \
                                        and  points(i,ii,iii)%r > rm
        # enddo
      
      else:
      # if the line of sight points to the moon's disk or the pixels with a 0-index
        forall(iii = -ni:ni)
          points(i,ii,iii)%compute =  False 
        # endforall
      # endif
    # enddo
  # enddo
      
  
# # end def line_of_sight
      
    
    
# dens is a 3d array, 3rd dimenssion is a def given
# at the uniformal grid (distance between the grid nodes = sampdist)
# with some missing values denoted by negative numbers
# the def integrates over the 3rd dimenssion using Euler's folmula
def Integral_over_LoS(dens0, image):
  use const
  implicit none
  #real, intent(in) :: dens0(-Hpix:Hpix, -Vpix:Vpix, -ni:ni)
  #real dens(-Hpix:Hpix, -Vpix:Vpix, -ni:ni)
  #real, intent(out) :: image(-Hpix:Hpix, -Vpix:Vpix)
  #integer i, ii, k
  #real tmpTr

  image = 00
  dens = dens0
  # double loop
  # over all the pixels
  # exclude the cross in the middle of the "image"
  for i  in range(-Hpix, Hpix):
  for ii  in range(-Vpix, Vpix):
    if(i  !=  0  and  ii  !=  0) :
      tmpTr = 00
      for k  in range(-ni, ni-1):
        tmpTr = tmpTr + dens(i,ii,k) + dens(i,ii,k+1)
      # enddo
      image(i,ii) = sampdist_large * 0.50 * tmpTr
    # endif
  # enddo
  # enddo  		

# # end def Integral_over_LoS  	



  
def result_image_out(image, m):
  use const
  implicit none
  #integer, intent(in) :: m
  #real, intent(in) :: image(-Hpix:Hpix,-Vpix:Vpix)
  #integer, parameter :: outchannel = 113
  #integer i 
  character(len = 55) fname
  character(len = 20) outformat
  character(len = 12) namesetup
  character(len = 2) tmp1, tmp2
  character(len = 1) tmp0
  
  write(tmp0,'(I1)') m
  
  tmp1 = 'I3'; tmp2 = 'I3'
  if(Vpix*2 < 100) tmp1 = 'I2'
  if(Vpix*2 < 10) tmp1 = 'I1'
  if(Hpix*2 < 100) tmp2 = 'I2'
  if(Hpix*2 < 10) tmp2 = 'I1'
  
  write(outformat,'(A1,I3,A1,I3,A11)') '(', (2*Vpix), 'x',(2*Hpix), '(ES12.4E2))'
  fname = './results/' // tmp0 // '.dat'
  
  open(outchannel, file = fname, status = 'replace')
  
    for i  in range(-Vpix, -1):
      write(outchannel,outformat) image(-Hpix:-1,i), image(1:Hpix,i)
    # enddo
    for i  in range(0, Vpix):
      write(outchannel,outformat) image(-Hpix:-1,i), image(1:Hpix,i)
    # enddo
  
  close(outchannel)
  write(*,*) 'result is in the file', fname

# # end def result_image_out  



# end module image_construction
