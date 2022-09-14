# PyPlumes
This is a workspace for a plume particle model in Python. <br />
Translated and modified from [this model](https://github.com/Veyza/dudi) originally written in Fortran. 

Run the testing/enceladus_test.py file for an example output with real data from Cassini flyby. A model for dust number density on Enceladus will be created and a plot comparing that with the real HRD observation during the Cassini flyby will be shown. 

# Tutorial 

See testing/enceladus_test.py or testing/europa_test.py for detailed examples with explanations.

The general instructions to run this package are as follows:

## Specify the key parameters 

open the file const.py and input values according to your chosen object:

- moon_mass                  mass of the moon, kg

- rm                         radius of the moon, meters


- rho                        density of the dust particles' material 
                             in kg/m^3,(needed if one wants to compute
                             mass density or mass fluxes)

- flux                       if set .TRUE. then dust flux through the
                             surface parallel to the moon surface is computed
                             instead of density

- p                          0 -- number density is computed, 1 -- mean
                             radius, 2 -- cross section, 3 -- mass density

- rmin                       lower boundary for function Gu(rmin, rmax),
                             microns

- rmax                       upper boundary for function Gu(rmin, rmax),
                             microns

- GRN                        number of Gu(Rmin,Rmax) values that are
                             precalculated other values are obtained
                             by interpolation from precalculated values

- order_R                    order of Gauss-Legendre quadrature formula
                             used for integration over particle radius R
                             (possible values are 5, 10, 20, 30)

- order_v_el                 order of Gauss-Legendre quadrature formula
                             used for integration over velocity to obtain
                             the particles density separately for particles
                             on bound (elliptic) trajectories
                             (possible values are 5, 10, 20, 30, 40, 50)
                        
- order_v_hy                 order of Gauss-Legendre quadrature formula
                             used for integration over velocity to obtain
                             the particles density separately for particles
                             on unbound (hyperbolic) trajectories
                             (possible values are 5, 10, 20, 30, 40, 50)       

##                        


