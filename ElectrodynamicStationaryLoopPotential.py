import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.integrate import simps
from scipy.integrate import quad
import copy
from matplotlib.colors import Normalize
import matplotlib.cm as cm

#########
##TO DO##
#########

#still fix the tuning of the static integral (i was too lazy to think about the parameterization)

#incorporate retarded time into the vector potential integral
##so this will be easy in a static current case
##i don't technically have to do it until i implement the loop rotation
##so this is basically done for this implementation

#implement the scalar potential (should I plot this as vector color? red to blue?)
##the charge has a domain of the loop coordinates, but the codomain is a constant scalar
##this should be an easier version of the vector potential since the loop is stationary

###########
##SUMMARY##
###########

#this uses the theory required for time evolution and electric-magnetic interaction

#there are now two integrals: one for electric scalar potential and one for magnetic vector potential
#the integrals at a time t rely on the current and charge distribution at a previous time t_r


##################
##IMPLEMENTATION##
##################

cmap = 'hot'

samples = 100

#initializing loop incline matrix
theta = np.radians(1)
c, s = np.cos(theta), np.sin(theta)
R_y = np.array(((c, 0.0, s),(0.0, 1.0, 0.0),(-s, 0.0, c)))

#gives you the coordinates of the stationary loop inclined at some angle (lamda in degrees)
#k is radius
def get_loop(lamda, k = 1.0):
    lamda = np.radians(lamda)
    return np.array(((k*np.cos(lamda)), (k*np.sin(lamda)), (0)))
#gives the current vector of the stationary loop inclined at some angle
#current is current magnitude
def get_current(lamda, current = 0.7):
    lamda = lamda+90
    lamda = np.radians(lamda)
    return np.array(((current*np.cos(lamda)), (current*np.sin(lamda)), (0)))


#gives you the coordinates along the rotated loop
#this is basically r(t)
def get_loop_rotated(t=0, steps = 100, k = 1.0):
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
        sol = R_z.dot(R_y.dot(get_loop(lamda, k=k)))
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

def get_vec_potential_at_r(r, current_dom, current_codom):
    to_be_integrated = np.zeros((current_codom.shape))
    for x in range(current_codom.shape[0]):
        #janky catch for division by zero
        cur_dist = get_dist(current_dom[x], r)
        if cur_dist != 0:
            to_be_integrated[x] = (current_codom[x][:])/cur_dist
        else:
            to_be_integrated[x] = np.array([0,0,0])

    return simps(to_be_integrated, dx=(3.14159*np.linalg.norm(current_dom[0])/samples), axis=0)

def get_scalar_potential_at_r(r, dom, codom):
    to_be_integrated = np.zeros((codom.shape))
    for x in range(codom.shape[0]):
        #janky catch for division by zero
        cur_dist = get_dist(dom[x], r)
        if cur_dist != 0:
            to_be_integrated[x] = (codom[x])/cur_dist
        else:
            to_be_integrated[x] = np.array([0,0,0])

    return simps(to_be_integrated, dx=(3.14159*np.linalg.norm(current_dom[0])/samples), axis=0)


#base data about loop coordinates and current vectors
x,y,z,current_dom = get_loop_rotated(steps=samples)
a,b,c,current_codom = get_current_rotated(steps=samples)

#boring charge map
charge = np.full((samples), 1.0)

base_x,base_y,base_z = ([] for i in range(3))
vec_x, vec_y, vec_z = ([] for i in range(3))

#i'm pretty sure that plugging this array in as colors just details the color of each vector in order as a scalar
scalar_pot = []

#test case
for cur_x in np.arange(-2.0, 2.0, .9):
    for cur_y in np.arange(-2.0, 2.0, .9):
        for cur_z in np.arange(-2.0, 2.0, .9):
        #for cur_z in [0.0]:
            cur_vec = get_vec_potential_at_r(np.array([cur_x,cur_y,cur_z]), current_dom, current_codom)
            print(cur_vec)

            scalar_pot.append(get_scalar_potential_at_r(np.array([cur_x,cur_y,cur_z]), current_dom, charge))
            
            base_x.append(cur_x)
            base_y.append(cur_y)
            base_z.append(cur_z)
            vec_x.append(cur_vec[0])
            vec_y.append(cur_vec[1])
            vec_z.append(cur_vec[2])


#this turns our scalar potential into valid color scalars
#norm = Normalize()
#norm.autoscale(scalar_pot) 
#this selects the color map for the scalars
colormap = cm.inferno
c = np.array(scalar_pot)
c = (c-c.min())/c.ptp()
c = np.concatenate((c, np.repeat(c,2)))
mask = charge==0
repeated_mask = np.concatenate((mask.ravel(), np.repeat(mask.ravel(),2)))
c = getattr(plt.cm, cmap)(c)
                              
fig = plt.figure()
ax = fig.gca(projection='3d')
#plotting the wire itself
ax.plot(x,y,z)
                              
#plotting the potential field
#not totally sure how colormap(norm(scalar_pot)) works
#p = ax.quiver(base_x, base_y, base_z, vec_x, vec_y, vec_z, cmap=colormap(norm(scalar_pot)))

p = ax.quiver(base_x, base_y, base_z, vec_x, vec_y, vec_z, cmap=cmap)

#ax.xaxis.set_ticks([])
#ax.yaxis.set_ticks([])
ax.set_aspect('equal')

p.set_array(np.linspace(np.min(charge), np.max(charge), 10))
fig.colorbar(p)
p.set_edgecolor(c)
p.set_facecolor(c)

plt.show()
