import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import variables as var
import integrator
import input
import output
import distributions as distf
import gu

r1 = 0.20
r2 = 20.0
europa_orbital_period = 3.06822*(10**5)	## in seconds
# number of sources
Ns = 4
#number of points on the SC trajectory for which the number density is to be calculated
nt = 80
	#integer i_s, i
	#real massflux(nt,2)
tnow = 0.0
	#real(8) dphi(nt), ddphi, tstep
	#real(8) mass_shallow, mass_steep, m1, m2
	#type(source_properties) source(Ns)
	#type(position_in_space) point(nt)
	
	call get_europa_input(Ns, source, nt, point, dphi)

	call mass_production(m1, source(1)%sd, r1, r2)
	call mass_production(m2, source(2)%sd, r1, r2)
	! integrate the production rate over time

	mass_steep = source(1)%production_rate * m2
	mass_shallow = source(1)%production_rate * m1
	
	write(*,*) '   '
	write(*,'(A84,e10.3,x,A2)') 'with the shallow size distribution the total &
								mass produced in a second', mass_shallow, 'kg'
	write(*,'(A84,e10.3,x,A2)') 'with the steep size distribution the total &
								mass produced in a second  ', mass_steep, 'kg'
	write(*,*) '   '
								
	do i_s = 1, Ns
		!$OMP PARALLEL PRIVATE(i) &
		!$OMP SHARED(i_s, point, source, massflux)
		!$OMP DO
		do i = 1, nt
			call DUDI(massflux(i,:), point(i), source(i_s), tnow)
		enddo
		!$OMP END DO
		!$OMP END PARALLEL
		call surface_deposition_out(i_s, massflux(:,1), nt, dphi)
	enddo

	
end
