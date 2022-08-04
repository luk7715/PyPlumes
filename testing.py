import numpy as np
import help 
import const as c 
import variables as var
import twobody as tb
import distributions as distf
import integrator as inte


t1 = np.array([2,4,6])
M1 = np.array([1,1,1])
A = np.zeros((3,3))
B = np.array([0,0,0])


A[0] = t1 
A[1] = ([1,3,10])
A[2] = ([2,1,1])
B[1] = np.dot(t1, M1)

Ainv = help.invert_matrix3(A)
#print (Ainv)

####need to test about file locations 
"""f = open("theta_geometry_output.txt","a")
f.write("\n")
f.write("\nTheta was found with an insufficient accuracy of"  + "from the geometry of hyperbola \n")
f.write("for r = \n" )
f.write("dphi = \n" )
f.close() """

class ejection_speed_properties:
    def __init__(self, ud_shape, umax, umin):
     self.ud_shape = ud_shape                 # parameter defining which distribution is used for ejection speed
     self.umax = umax                         # gas velocity
     self.umin = umin                         # a parameter in velocity distribution
    
    def get_ud_shape(self):
     return self.ud_shape             

test_ud = ejection_speed_properties(10,100,1)
print(test_ud.ud_shape)
print(test_ud.umax)
print(test_ud.umin)


class source_properties(ejection_speed_properties):                     # parameters of dust ejection
    def __init__(self, ejection_speed_properties):
#ype(ejection_speed_properties) ud       # parameters of ejection speed distribution
        self.ud_shape = ejection_speed_properties.ud_shape
        self.ud_umax = ejection_speed_properties.umax
        self.ud_umin = ejection_speed_properties.umin

source = source_properties(test_ud)


test_source = source
print (source.ud_umax)

