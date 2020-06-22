import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.integrate import simps
import copy

#########
##TO DO##
#########

#iron out the bug causing incorrect vector field
#should loosely follow direction of loop (right?)
#magnitudes look approximately right...
#the issue seems to be in the integration. are we not integrating over the whole loop?

###########
##SUMMARY##
###########

#The magnetostatic case allows us to calculate the magnetic potential by
#taking an integral over the current.

#The integral is something like A(r) = (mu_naught/(4*pi))*Integral[I(cur)/dist(r,cur) d(cur)] where I is the current vector at the point and cur is an infinitesmal of the curve we're integrating over.

#I'm gonna try and keep this general by using some kind of integration method and sampling many r points

##################
##IMPLEMENTATION##
##################

#radius of loop
k = 1.0
#affects current magnitude (Not exactly current magnitude; i'm lazy)
current = .1

#initializing loop incline matrix
theta = np.radians(1)
c, s = np.cos(theta), np.sin(theta)
R_y = np.array(((c, 0.0, s),(0.0, 1.0, 0.0),(-s, 0.0, c)))

#gives you the coordinates of the stationary loop inclined at some angle (lamda in degrees)
def get_loop(lamda):
    lamda = np.radians(lamda)
    return np.array(((k*np.cos(lamda)), (k*np.sin(lamda)), (0)))
#gives the current vector of the stationary loop inclined at some angle
def get_current(lamda):
    lamda = lamda+90
    lamda = np.radians(lamda)
    return np.array(((current*np.cos(lamda)), (current*np.sin(lamda)), (0)))


#gives you the coordinates along the rotated loop
#this is basically r(t)
def get_loop_rotated(t=0, steps = 100):
    #getting the rotation matrix for time t
    t = np.radians(t)
    c_1, s_1 = np.cos(t), np.sin(t)
    R_z = np.array(((c_1, -s_1, 0.0),(s_1, c_1, 0.0),(0.0, 0.0, 1.0)))

    x = []
    y = []
    z = []
    vectors = []
    #getting the bits of the loop at the desired resolution
    step_size = 360.0/steps
    for cur in range(steps):
        lamda = step_size*cur
        sol = R_z.dot(R_y.dot(get_loop(lamda)))
        x.append(sol[0])
        y.append(sol[1])
        z.append(sol[2])
        vectors.append(np.array(sol))
        
    return x,y,z,np.stack(vectors, axis=0)

#gives you the vectors along the rotated loop
#this is 
def get_current_rotated(t=0, steps = 100):
    #getting the rotation matrix for time t
    t = np.radians(t)
    c_1, s_1 = np.cos(t), np.sin(t)
    R_z = np.array(((c_1, -s_1, 0.0),(s_1, c_1, 0.0),(0.0, 0.0, 1.0)))

    x = []
    y = []
    z = []
    vectors = []
    #getting the bits of the loop at the desired resolution
    step_size = 360.0/steps
    for cur in range(steps):
        lamda = step_size*cur
        sol = R_z.dot(R_y.dot(get_current(lamda)))
        x.append(sol[0])
        y.append(sol[1])
        z.append(sol[2])
        vectors.append(np.array(sol))
        
    return x,y,z,np.stack(vectors, axis=0)

#vectorized function for crunching the inside of the integral, J/|r-cur|
def get_dist(current_dom, r):
    return np.linalg.norm(current_dom-r)

def get_potential_at_r(r, current_dom, current_codom):
    to_be_integrated = np.zeros((current_codom.shape))
    for x in range(current_codom.shape[0]):
        #janky catch for division by zero
        cur_dist = get_dist(current_dom[x], r)
        if cur_dist != 0:
            to_be_integrated[x] = (current_codom[x][:])/cur_dist
        else:
            to_be_integrated[x] = np.array([0,0,0])

    print("TO BE INTEGRATED")
    print(to_be_integrated)
    print("CURRENT DOMAIN")
    print(current_dom)
    return simps(to_be_integrated, x=current_dom, axis=0)

def debug_to_be_integrated(r, current_dom, current_codom):
    to_be_integrated = np.zeros((current_codom.shape))
    for x in range(current_codom.shape[0]):
        #janky catch for division by zero
        cur_dist = get_dist(current_dom[x], r)
        if cur_dist != 0:
            to_be_integrated[x] = (current_codom[x][:])/cur_dist
        else:
            to_be_integrated[x] = np.array([0,0,0])
    #print(to_be_integrated)
    return to_be_integrated

#base data about loop coordinates and current vectors
x,y,z,current_dom = get_loop_rotated(steps=100)
a,b,c,current_codom = get_current_rotated(steps=100)

base_x,base_y,base_z = ([] for i in range(3))
vec_x, vec_y, vec_z = ([] for i in range(3))

#print("Current domain values")
#print(current_dom)
#print("Current codomain values")
#print(current_codom)

#test case
for cur_x in np.arange(-2.0, 2.0, .7):
    for cur_y in np.arange(-2.0, 2.0, .7):
        for cur_z in np.arange(-2.0, 2.0, .7):
        #for cur_z in [0.0]:
            cur_vec = get_potential_at_r(np.array([cur_x,cur_y,cur_z]), current_dom, current_codom)
            print(cur_vec)
            base_x.append(cur_x)
            base_y.append(cur_y)
            base_z.append(cur_z)
            vec_x.append(cur_vec[0])
            vec_y.append(cur_vec[1])
            vec_z.append(cur_vec[2])

            
fig = plt.figure()
ax = fig.gca(projection='3d')
#plotting the wire itself
ax.plot(x,y,z)

#plotting the current vectors
#this looks right
#ax.quiver(x, y, z, a, b, c)

#debug test
#test_data = debug_to_be_integrated(np.array([1.5,0.0,0.0]), current_dom, current_codom)
#ax.quiver(x, y, z, test_data[:,0], test_data[:,1], test_data[:,2])

#plotting the potential field
#this does not look right
ax.quiver(base_x, base_y, base_z, vec_x, vec_y, vec_z)

plt.show()
