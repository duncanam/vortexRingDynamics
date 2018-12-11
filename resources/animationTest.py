import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3


######################################
# GIVEN INFO, CONSTANTS AND OTHER VARIABLES

pi = np.pi #pi
R0 = 1 #initial radius of both rings
H0 = 1 #initial spacing of the rings
Gamma = 1 #vortex strength for both rings
res = 100 #nubmer of points used to model each vortex ring
global R
R = R0
######################################

theta = np.linspace(0,2*pi,res) #defining a vector of angles to make a circle
y1 =  R0*np.cos(theta) #x values of points for the first ring (BACK)
z1 =  R0*np.sin(theta) #y values of points for the first ring (BACK)


def animate(i):
    global R
    R = R-0.01
    if(R<0.01):
        R = R0
    line.set_xdata(R*np.cos(theta))  # update the data.
    line.set_ydata(R*np.sin(theta))



    return line


fig, ax = plt.subplots()
line, = ax.plot(y1,z1)



ani = animation.FuncAnimation(
    fig, animate, interval=2, blit=False, save_count=50)

# To save the animation, use e.g.
#
# ani.save("movie.mp4")
#
# or
#
# from matplotlib.animation import FFMpegWriter
# writer = FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save("movie.mp4", writer=writer)
plt.axis("equal")
plt.show()
