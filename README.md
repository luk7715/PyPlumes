# PyPlumes
This is a workspace for a plume particle model in Python. <br />
Translated and modified from [this model](https://github.com/Veyza/dudi) originally written in Fortran. 

Run the testing/enceladus_test.py file for an example output with real data from Cassini flyby. A model for dust number density on Enceladus will be created and a plot comparing that with the real HRD observation during the Cassini flyby will be shown. 

# Tutorial 

See testing/enceladus_test.py or testing/europa_test.py for detailed examples with explanations.

The general instructions to run this package are as follows:

## Specify the key parameters 
### ALWAYS CHECK THE const.py FILE TO HAVE CORRECT VALUES 

open the file const.py and input values according to your chosen object:

- **moon_mass**: mass of the moon, kg

- **rm**: radius of the moon, meters


- **rho**: density of the dust particles' material 
                             in kg/m^3,(needed if one wants to compute
                             mass density or mass fluxes)

- **flux**: if set .TRUE. then dust flux through the
                             surface parallel to the moon surface is computed
                             instead of density

- **p**: 0 -- number density is computed, 1 -- mean
                             radius, 2 -- cross section, 3 -- mass density

- **rmin**: lower boundary for function Gu(rmin, rmax),
                             microns

- **rmax**: upper boundary for function Gu(rmin, rmax),
                             microns

- **GRN**: number of Gu(Rmin,Rmax) values that are
                             precalculated other values are obtained
                             by interpolation from precalculated values

- **order_R**: order of Gauss-Legendre quadrature formula
                             used for integration over particle radius R
                             (possible values are 5, 10, 20, 30)

- **order_v_el**: order of Gauss-Legendre quadrature formula
                             used for integration over velocity to obtain
                             the particles density separately for particles
                             on bound (elliptic) trajectories
                             (possible values are 5, 10, 20, 30, 40, 50)
                        
- **order_v_hy**: order of Gauss-Legendre quadrature formula
                             used for integration over velocity to obtain
                             the particles density separately for particles
                             on unbound (hyperbolic) trajectories
                             (possible values are 5, 10, 20, 30, 40, 50)       


##  Specify the ejection distribution functions
Go to file "distributions.py" for functions describing the size, ejection speed,  ejection direction distributions, and the time-dependent dust production rate function. If necessary, users can write their desired distribution in the corresponding functions. "Match" operator is used to select the most suitable case.

Detailed documentations and instructions are in file "distributions.py"

##  Supply input data
All the input data are saved in file "variables.py" into their corresponding classes (source, point, and ejection distribution) so that they can be easily accessed and refrenced by other functions. 
See "variables.py" for detailed description of each class and variable.

To read in data, go to file "input.py". There are existing example functions for constructing the aforementioned classes. When necessary, the users can write their own data input functions in this file. See file "input.py for detailed documentations. 

## Call Function DUDI 
Function DUDI(short for Dust Distribution) in the file "integrator.py" is the code that performs numerical integrations to compute density at the given point with known source parameters at the time tnow. It is the core of this model. 

The file "template.py" is provided as a template of the main program. It calls the subroutines from the module input_data to obtain an array of sources ejecting dust and an array of the points in space where one wants to compute density. Then the program makes a double loop over the sources and over the points calling the subroutine DUDI from the module integrator. In each point the result density is the sum of densities of the dust from all the sources. OpenMP is used to speed up the calculations as for each pair of source and point the calculations are independent. Finally, the main program calls the function from module dataout to write the result into the file.

Modify the file "template.py" and the subroutines for input and output for your own needs.

