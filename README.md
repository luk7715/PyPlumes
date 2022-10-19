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

## Call main working function
Function DUDI(short for Dust Distribution) in the file "integrator.py" is the code that performs numerical integrations to compute density at the given point with known source parameters at the time tnow. It is the core of this model. 
Make sure to cover all sources and desired points with DUDI.

## Output of the results
The function result_out in the module "output.py" is provided for outputing the results. It writes the result into the file "twobody_model_result.txt" in the directory "./results/" as follows:

- 1st column is the density at the point

- 2nd column is the density of particles on elliptic orbits

- 3rd column is the density of particles on hyperbolic orbits
	
- 4th column is the spacecraft's radial distance in meters

- 5th column is the spacecraft's latitude in degrees

- 6th column is the spacecraft's eastern longitude in degrees

1st column is sum of the 2nd and the 3rd.

The users can write their own output functions. 



# Examples and Templates

## Template
The file "template.py" is provided as a template of the main program. It calls the functions from the module input to obtain an array of sources ejecting dust and an array of the points in space where one wants to compute density. Then the program makes a double loop over the sources and over the points calling the function DUDI from the module integrator. In each point the result density is the sum of densities of the dust from all the sources. Finally, the main program calls the function from module output to write the result into the file.

Modify the file "template.py" and the functions for input and output for your own needs.

## Example 1. The number density profile of the E2 flyby of the Cassini spacecraft at Enceladus

The file "./testing/enceladus_test.py" performs a loop over 100 points along the trajectory of Cassini and computes in these points the number density of dust from one tilted jet representing the Enceladus dust plume.

Set the parameters in const.py to the following values:

moon_mass = 1.08022*(10**20)       #Enceladus's mass in kg 

rm = 252*(10**3)                   #Enceladus's mean radius in meters

flux = .FALSE.                     #we are currently not interested in flux

p = 0                              #this will give us a number density profile 

rmin = 1.60                        #micron, the lower threshold of the HRD sensitivity 

rmax = 6.00                        #micron, a fairly large upper boundary, it is assumed that the probability to detect a particle larger than rmax is negligibly small

GRN = 2000                         #number of precalculated values of Gu integral 

order_R = 30                       #order of integration Gu(Rmin,Rmax) over R to obtain Gu 

order_v_el = 50                    #order of Gaussian quadrature for integration of n(r, alpha, beta, v, theta, lambda) over velocity interval corresponding to elliptic trajectories 

order_v_hy = 20                    #order of Gaussian quadrature for integration of n(r, alpha, beta, v, theta, lambda) over velocity interval corresponding to hyperbolic trajectories

The directory input_data_files contains files "Enceladus_jet.dat" and "Cassini_E2_flyby.dat" with the input data. "Enceladus_jet.dat" contains a line with the source properties. The coordinates, zenith angles and azimuth of the jet were adopted from the paper of Porco et al., 2014.

"Cassini_E2_flyby.dat" contains the coordinates of the Cassini spacecraft at the required moments. The coordinates were obtained with the SPICE package and stored into this file for simplicity.


Run the script in "./testing/enceladus_test.py" will construct the model and output results into the file "./results/E2_profile.dat". The output consists of 2 columns. The first column is the time in seconds counted from 2005-07-14 19:49:21 (the moment of the closest approach), the second column is the dust number-density. A plot with comparison of the model profile with the observational data will be generated and shown on screen

## Example 2. Dust deposition on the surface of Europa

The file "./testing/europa_test.py" performs calculation of the dust flux onto the surface of Europa. The total mass production rate is computed for 2 different size distributions. The surface deposition is calculated separately for 4 sources different in size and ejection direction distributions

Set the parameters in const.py to the following values:

moon_mass = 4.8*(10*22)             #Europa's mass in kg 

rm = 3.1216*(10**6)                 #Europa's mean radius in meters

rho = 920.0                         #dust grains density in kg/m^3  (assuming pure water ice composition) 

flux = .TRUE.                       #this time we are calculating flux

p = 3                               #this will give us a mass flux 

rmin = 0.20                         #the lower limit of the applied size distributions, microns 

rmax = 20.0                         #the upper limit of the applied size distributions, microns

GRN = 2000                          #number of precalculated values of Gu integral 

order_R = 10                        #order of integration Gu(Rmin,Rmax) over R to obtain Gu 

order_v_el = 30                     #order of Gaussian quadrature for integration of n(r, alpha, beta, v, theta, lambda) over velocity  interval corresponding to elliptic trajectories

order_v_hy = 5                      #does not matter in this example because all the particles are on elliptic trajectories

No input files are used. The parameters of the sources and the points on the surface, where the flux is calculated, are set in the function get_europa_input in the file "input.py


Run the file "./testing/europa_test.py" to produce a model for europa dust deposition. The result will be put in the files "narrow_jet_shallow_sd.dat", "diffuse_source_steep_sd.dat", "diffuse_source_shallow_sd.dat", and "narrow_jet_steep_sd.dat" in the directory "./results". They correspond to the four sources with different combinations of the size distribution and the ejection direction distribution. Each output file consists of two columns where the first column will be the distance from the source in km, and the second column will be the mass flux in kg/m^2/s. The total mass production rate for each size distribution will be printed in the terminal. 
A plot will be produced and saved to the file "mass_deposition.png" in the directory "./results".