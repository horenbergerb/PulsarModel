import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.integrate import simps

#########
##TO DO##
#########

#Crunch what the line integral actually looks like
#Parameterize the curve and whatnot
#Don't get frustrated

###########
##SUMMARY##
###########

#The magnetostatic case allows us to calculate the magnetic potential by
#taking a line integral over the current.

#The integral is something like A(r) = (mu_naught/(4*pi))*Integral[I(cur)/dist(r,cur) d(cur)] where I is the current vector at the point and cur is an infinitesmal of the curve we're integrating over.

#I'm gonna try and keep this general by using some kind of integration method and sampling many r points

##################
##IMPLEMENTATION##
##################

#radius of loop
k = 3.0

#end time of simulation (in degrees; don't question this)
t_f = 360

#initializing loop incline matrix
theta = np.radians(30)
c, s = np.cos(theta), np.sin(theta)
R_y = np.array(((c, 0.0, s),(0.0, 1.0, 0.0),(-s, 0.0, c)))

#gives you the coordinates of the loop at some angle (lamda in degrees)
def get_loop(lamda):
    lamda = np.radians(lamda)
    return np.array(((k*np.cos(lamda)), (k*np.sin(lamda)), (0)))

def get_current(lamda):
    lamda = lamda+90
    lamda = np.radians(lamda)
    return np.array(((k*np.cos(lamda)), (k*np.sin(lamda)), (0)))


#gives you the coordinates along the rotated loop
#this is basically r(t)
def get_loop_rotated(t, steps = 100):
    #getting the rotation matrix for time t
    t = np.radians(t)
    c, s = np.cos(t), np.sin(t)
    R_z = np.array(((c, -s, 0.0),(s, c, 0.0),(0.0, 0.0, 1.0)))

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
        
    return x,y,z,vectors

#gives you the vectors along the rotated loop
#this is 
def get_current_rotated(t, steps = 100):
    #getting the rotation matrix for time t
    t = np.radians(t)
    c, s = np.cos(t), np.sin(t)
    R_z = np.array(((c, -s, 0.0),(s, c, 0.0),(0.0, 0.0, 1.0)))

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
        
    return x,y,z,vectors

#vectorized function for crunching the inside of the integral, J/|r|
def divide_by_distance(r, current_dom, current_codom):
    return current_codom/np.linalg.norm(current_dom-r)

def get_potential_at_r(r, current_dom, current_codom):
    target_func = np.vectorize(divide_by_distance)
    to_be_integrated = target_func(r, current_dom, current_codom)


    return simps(to_be_integrated, x=current_dom, axis=0)
    

#initializing for the animation
x,y,z,current_dom = get_loop_rotated(0, steps=100)
a,b,c,current_codom = get_current_rotated(0, steps=100)

get_potential_at_r(np.array([1,1,1]), current_dom, current_codom)

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#line, = ax.plot(x,y,z)
#ax.plot(x,y,zs=z)
