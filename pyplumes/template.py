# This file is a part of PyPlumes, a Python model for dust ejection dynamics
# on atmosphereless surfaces
# Version 1.0.0
# Translated and modified from the Fortran-95 code in the following paper
# Anastasiia Ershova and Jürgen Schmidt, 
# Two-body model for the spatial distribution of dust ejected from
# an atmosphereless body, 2021, A&A, 650, A186 
# File: template.py
# Description: A template for applying this model

# Author: Eulrika(Yiqi) Wu
# E-mail: ulkw517@g.ucla.com

import numpy as np
import matplotlib.pyplot as plt
import variables as var
import integrator
import input
import output


Ns = 1                                     # number of sources
nt = 1                                     # number of points on the SC trajectory for which the number density is to be calculated
tnow = 600.0                               #time for the calculation
density = np.zeros([nt,2])
tmp_res = np.zeros([nt,2])

##Input data
var.source = input.read_sources_params('./input_data_files/sources.dat', Ns)
var.point = input.read_spacecraft_coordinates('./input_data_files/points.dat', nt)


#loop through all sources
for i in range (0, Ns):
    for x in range (0,nt): 
        #loop through all data points   
        var.point.r = var.point.r[x]; 
        var.point.alpha = var.point.alpha[x]
        var.point.beta = var.point.beta[x]; var.point.r_scaled = var.point.r_scaled[x]
        var.point.rvector = var.point.rvector[x,:]
        
        #run the integration functions for each point to get the corresponding density 
        tmp_res[x,:] = integrator.DUDI(tnow)
    density = density + tmp_res

##Output resulat to files
output.result_out(density,nt)


