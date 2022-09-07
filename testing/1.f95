! the subroutine provides the input data for the example
! calculations of the Europa surface depositions
subroutine get_europa_input(Ns, sources, nt, points, dphi)
	use define_types
	use gu
	implicit none
	integer, intent(in) :: Ns, nt
	real(8), intent(out) :: dphi(nt)
	type(source_properties), intent(out) :: sources(Ns)
	
	type(position_in_space), intent(out) :: points(nt)
	integer i, ead_choice(4), sd_choice(4)
	
	ead_choice = (/1, 3, 3, 1/)
	sd_choice = (/2, 3, 2, 3/)
	
	! define 4 sources with the same coordinates and verticle axis
	! of symmetry but different size- and ejection direction distributions
	do i = 1, Ns
		sources(i)%alphaM = halfpi
		sources(i)%betaM = 0d0
		sources(i)%zeta = 0d0
		sources(i)%eta = 0d0
		sources(i)%production_fun = 0
		sources(i)%production_rate = 1d14
		sources(i)%ud%ud_shape = 1
		sources(i)%ud%umin = 0d0
		sources(i)%ud%umax = 500d0
		sources(i)%ejection_angle_distr = ead_choice(i)
		sources(i)%sd = sd_choice(i)
		sources(i)%r = rm

						
		sources(i)%rrM(1) = rm * sin(sources(i)%alphaM) &
								* cos(sources(i)%betaM)
		sources(i)%rrM(2) = rm * sin(sources(i)%alphaM) &
								* sin(sources(i)%betaM)	
		sources(i)%rrM(3) = rm * cos(sources(i)%alphaM)
			
		sources(i)%is_jet = .TRUE.
		
		call jet_direction(sources(i)%betaM, sources(i)%zeta, sources(i)%eta, &
							sources(i)%rrM, sources(i)%symmetry_axis)
		call Gu_integral(sources(i)%ui, sources(i)%Gu_precalc, sources(i)%sd, &
							sources(i)%ud)
	enddo
	
	! the points are equidistantly placed on a 10deg arc
	! having the source at one of the arc's ends
	forall(i = 1:nt) dphi(i) = 2d0 * dble(i) / dble(nt) * deg2rad
	
	do i = 1, nt
		points(i)%r = rm
		points(i)%alpha = halfpi 
		points(i)%beta = dphi(i)
		points(i)%rvector(1) = points(i)%r * sin(points(i)%alpha) &
											* cos(points(i)%beta)
		points(i)%rvector(2) = points(i)%r * sin(points(i)%alpha) &
											* sin(points(i)%beta)
		points(i)%rvector(3) = points(i)%r * cos(points(i)%alpha)
		points(i)%r_scaled = 1d0
		points(i)%compute = .TRUE.
	enddo				


end subroutine get_europa_input

