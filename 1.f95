! input the parameters of the sources from the file with the name
! stored in the variable fname
subroutine read_sources_params(fname,)
	use const
	use define_types
	use gu
	implicit none
	integer, intent(in) :: Ns
	type(source_properties), intent(out) :: sources(Ns)
	integer i
	character(*), intent(in) :: fname
	
	open(100, file = fname, status = 'old')
		do i = 1, Ns
			read(100,*) sources(i)%alphaM, sources(i)%betaM, &
						sources(i)%zeta, sources(i)%eta, &
						sources(i)%production_fun, sources(i)%production_rate, &
						sources(i)%ud%ud_shape, sources(i)%ud%umin, &
						sources(i)%ud%umax, &
						sources(i)%ejection_angle_distr, &
						sources(i)%sd
			sources(i)%r = rm
			sources(i)%alphaM = sources(i)%alphaM * deg2rad
			sources(i)%betaM = sources(i)%betaM * deg2rad
			sources(i)%alphaM = halfpi - sources(i)%alphaM
			
			sources(i)%zeta = sources(i)%zeta * deg2rad
			sources(i)%eta = sources(i)%eta * deg2rad
						
			sources(i)%rrM(1) = rm * sin(sources(i)%alphaM) &
									* cos(sources(i)%betaM)
			sources(i)%rrM(2) = rm * sin(sources(i)%alphaM) &
									* sin(sources(i)%betaM)
			sources(i)%rrM(3) = rm * cos(sources(i)%alphaM)
			
			sources(i)%is_jet = .TRUE.
			
			call jet_direction(sources(i)%betaM, sources(i)%zeta, &
					sources(i)%eta, sources(i)%rrM, sources(i)%symmetry_axis)
					
			call Gu_integral(sources(i)%ui, sources(i)%Gu_precalc, &
							sources(i)%sd, sources(i)%ud)
									
		enddo
	close(100)

end subroutine read_sources_params