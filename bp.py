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
gamma = np.array([1,1,2,1.5,1,0.5,0.8,1,1.2,5,5,2,2,2,2,2,2,2,1.6,1,1,1,1,0.5])
radius = np.array([1,2,1,1.5,1,1,2,0.8,0.7,2.0,2,4,5,6,7,7,7.2,6.8,6.5,3,3,1.6,2.5,0.8]) #np.linspace(0.1,10,1000)
height = np.array([1,1,1,1,0.5,0.5,0.2,0.4,0.6,1.0,2,1,1,1,1,1.2,1.4,1.4,1.5,0.2,1,0.2,1,1]) #np.linspace(10,20,1000)
T = np.array([12.81,17.6,6.61,10.81,4.61,8.82,1.21,3.01,4.41,3.62,10.41,10.01,10.02,10.21,10.21,14.4,19.62,19.61,27.6,0.841,19.01,0.82,18.82,21.62])

#####################################################
# CALCULATE THE PLOT GROUPS
p1 = pi1(gamma, T, radius, height)
p2 = pi2(radius, height)

#####################################################
# Curve fitting
def exponenial_func(x, a, b, c):
    return a*np.exp(-b*x)+c


popt, pcov = curve_fit(exponenial_func, p2, p1, p0=(1, 1e-6, 1))
xx = np.linspace(0, 15, 100)
yy = exponenial_func(xx, *popt)

#Print the exponential fit results:
print('pi1 = a*exp(-b*pi2) + c')
print(popt)
print('a =',popt[0])
print('b =',popt[1])
print('c =',popt[2])
#####################################################
# PLOT THE RESULTS

plt.figure()
plt.plot(p2,p1,'.',markersize="10")
plt.plot(xx,yy,'r')
plt.title('Leapfrogging Vortex Data Plotted with Curve Fit')
plt.xlabel('R/H')
plt.ylabel('GT/(RH)')
plt.show()
