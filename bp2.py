#####################################################
# ASEN 5051
# Project 2
# Buckingham Pi Plotting Script
# Lucas Calvert and Duncan McGough
#####################################################
# IMPORT PACKAGES
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#####################################################
# DEFINE DIMENSIONLESS GROUPS
def pi1(g,t,r,h):
    return g*t/(r*h)

def pi2(r,h):
    return r/h

#####################################################
# POPULATE VARIABLES
gamma = np.array([1,1,2,1.5,1,0.5,0.8,1,1.2,5,5,2,2,2,2,2,2,2,1.6,1,1,1,0.5])
radius = np.array([1,2,1,1.5,1,1,2,0.8,0.7,2.0,2,4,5,6,7,7,7.2,6.8,6.5,3,1.6,2.5,0.8]) #np.linspace(0.1,10,1000)
height = np.array([1,1,1,1,0.5,0.5,0.2,0.4,0.6,1.0,2,1,1,1,1,1.2,1.4,1.4,1.5,1,0.2,1,1]) #np.linspace(10,20,1000)
T = np.array([1.01,2.81,0.51,1.21,0.75,1.41,2.41,0.45,0.405,0.561,0.841,4.4,7.01,10.4,14.3,13.8,13.8,12.7,14.3,5.4,1.21,4.0,1.61])

#####################################################
# CALCULATE THE PLOT GROUPS
p1 = pi1(gamma, T, radius, height)
p2 = pi2(radius, height)

#####################################################
# Curve fitting
fit = np.polyfit(p2,p1,1)

##Print the linera fit results:
print('pi1 = a*p2 + b')
print(fit)
print('a =',fit[0])
print('b =',fit[1])

xFit = np.linspace(0,np.max(p2)+1,100)
yFit = fit[0]*xFit + fit[1]
#####################################################
# PLOT THE RESULTS

plt.figure()
plt.plot(p2,p1,'.',markersize="10")
plt.plot(xFit,yFit)
#plt.plot(xx,yy,'r')
plt.title('Leapfrogging Vortex Data Plotted with Curve Fit, H=H/2')
plt.xlabel('R/H')
plt.ylabel('GT/(RH)')
plt.show()
