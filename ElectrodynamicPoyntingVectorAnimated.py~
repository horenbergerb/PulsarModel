import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.integrate import simps
from scipy.integrate import quad
import copy
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from scipy.constants import c
from matplotlib import animation

#Color stuff was implemented using this:
#https://stackoverflow.com/questions/28420504/adding-colors-to-a-3d-quiver-plot-in-matplotlib
#the last answer was the only solution I could find that allowed for a colorbar

#########
##TO DO##
#########

#still fix the tuning of the static integral (i was too lazy to think about the parameterization)

#optimize calculations by reducing calls to get_current_rotated and get_loop_rotated

#implement animation that generates frames using the code from the pre-animated version

###########
##SUMMARY##
###########

#this uses the theory required for time evolution and electric-magnetic interaction

#there are now two integrals: one for electric scalar potential and one for magnetic vector potential
#the integrals at a time t rely on the current and charge distribution at a previous time t_r


##################
##IMPLEMENTATION##
##################

#graphics initialization

t_f = 360.0
frames = 20

quiver = None
cmap = None

fig = plt.figure()
ax = fig.gca(projection='3d')

cmap = 'hot'

#speed_light = 1.0
speed_light = c

#number of times to sample at within one revolution
time_samples = 100

samples = 100

#initializing loop incline matrix
theta = np.radians(45)
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
        
    return np.stack(vectors, axis=0)

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
        
    return np.stack(vectors, axis=0)

#vectorized function for crunching the inside of the integral, J/|r-cur|
def get_dist(current_dom, r):
    return np.linalg.norm(current_dom-r)

def get_vec_potential_at_r(r, t=0):
    to_be_integrated = None
    ###########################################
    #Get current via retarded time config here#
    ###########################################
    current_dom = get_loop_rotated(t=t)
    current_codom = get_current_rotated(t=t)
    to_be_integrated = np.zeros((current_codom.shape))

    for x in range(current_codom.shape[0]):
        cur_dist = get_dist(current_dom[x], r)
        retarded_current_dom = get_loop_rotated(t=t-(cur_dist/speed_light))
        retarded_current_codom = get_current_rotated(t=t-(cur_dist/speed_light))
        #janky catch for division by zero
        if cur_dist != 0:
            to_be_integrated[x] = (retarded_current_codom[x][:])/cur_dist
        else:
            to_be_integrated[x] = np.array([0,0,0])

    return simps(to_be_integrated, dx=(3.14159*np.linalg.norm(retarded_current_dom[0])/samples), axis=0)

#NOTE: Assumes constant charge along wire
def get_scalar_potential_at_r(r, charge, t=0):
    ###########################################
    #Get current via retarded time config here#
    ###########################################
    dom = get_loop_rotated(t=t)
    codom = charge
    to_be_integrated = np.zeros((codom.shape))
    for x in range(codom.shape[0]):
        cur_dist = get_dist(dom[x], r)
        retarded_dom = get_loop_rotated(t=t-(cur_dist/speed_light))
        retarded_codom = charge
        #janky catch for division by zero
        if cur_dist != 0:
            to_be_integrated[x] = (retarded_codom[x])/cur_dist
        else:
            to_be_integrated[x] = np.array([0,0,0])

    return simps(to_be_integrated, dx=(3.14159*np.linalg.norm(retarded_dom[0])/samples), axis=0)

#time sample domain
times = np.linspace(0, 360, 360/time_samples)

def plot_at_time(t, line):
    #have to recreate these for each frame
    global quiver
    global cmap
    global fig
    global ax
    
    #base data about loop coordinates and current vectors
    current_dom = get_loop_rotated(steps=samples, t=t)
    current_codom = get_current_rotated(steps=samples, t=t)

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
                cur_vec = get_vec_potential_at_r(np.array([cur_x,cur_y,cur_z]), t=t)
                #print(cur_vec)

                scalar_pot.append(get_scalar_potential_at_r(np.array([cur_x,cur_y,cur_z]), charge, t=t))
            
                base_x.append(cur_x)
                base_y.append(cur_y)
                base_z.append(cur_z)
                vec_x.append(cur_vec[0])
                vec_y.append(cur_vec[1])
                vec_z.append(cur_vec[2])

    colormap = cm.inferno
    c = np.array(scalar_pot)
    c = (c-c.min())/c.ptp()
    c = np.concatenate((c, np.repeat(c,2)))
    mask = charge==0
    repeated_mask = np.concatenate((mask.ravel(), np.repeat(mask.ravel(),2)))
    c = getattr(plt.cm, cmap)(c)

    #plotting the wire itself
    line.set_data(current_dom[:,0], current_dom[:,1])
    line.set_3d_properties(current_dom[:,2])
    
    #plotting the potential field
    if quiver != None:
        quiver.remove()
    quiver = ax.quiver(base_x, base_y, base_z, vec_x, vec_y, vec_z, cmap=cmap)

    ax.set_aspect('equal')

    quiver.set_array(np.linspace(np.min(charge), np.max(charge), 10))
    #only add colorbar in first iteration
    if t == 0:
        fig.colorbar(quiver)
    quiver.set_edgecolor(c)
    quiver.set_facecolor(c)
    print("Frame rendering complete")

#initialize the loop
current_dom = get_loop_rotated(steps=samples, t=0.0)
line, = ax.plot(current_dom[:,0], current_dom[:,1], current_dom[:,2])
#initialize the quiver
#plot_at_time(0.0, line,)

ani = animation.FuncAnimation(fig, plot_at_time, (x*(t_f/frames) for x in range(0, frames)), fargs=(line,), interval=200, blit=False)
ani.save('field.gif', writer='imagemagick')
plt.show()
