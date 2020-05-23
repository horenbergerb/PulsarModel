import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import animation

#radius of loop
k = 3.0

#end time of simulation (in degrees; don't question this)
t_f = 360

#initializing loop angle matrix
theta = np.radians(30)
c, s = np.cos(theta), np.sin(theta)
R_y = np.array(((c, 0.0, s),(0.0, 1.0, 0.0),(-s, 0.0, c)))

#initializing the loop itself
def get_loop(lamda):
    lamda = np.radians(lamda)
    return np.array(((k*np.cos(lamda)), (k*np.sin(lamda)), (0)))

def get_loop_rotated(t, steps = 100):
    #getting the rotation matrix for time t
    t = np.radians(t)
    c, s = np.cos(t), np.sin(t)
    R_z = np.array(((c, -s, 0.0),(s, c, 0.0),(0.0, 0.0, 1.0)))


    x = []
    y = []
    z = []
    #getting the bits of the loop at the desired resolution
    step_size = 360.0/steps
    for cur in range(steps):
        lamda = step_size*cur
        sol = R_z.dot(R_y.dot(get_loop(lamda)))
        x.append(sol[0])
        y.append(sol[1])
        z.append(sol[2])
        
    return x,y,z

#update function for the animator
def update(num, line):
    x,y,z = get_loop_rotated(num)
    line.set_data(x,y)
    line.set_3d_properties(z)

#initializing for the animation
x,y,z = get_loop_rotated(0, steps=100)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
line, = ax.plot(x,y,z)
#ax.plot(x,y,zs=z)

ani = animation.FuncAnimation(fig, update, t_f, fargs=(line,), interval = 1000/t_f, blit=False)
ani.save('spinningloop.gif', writer='imagemagick')
plt.show()
