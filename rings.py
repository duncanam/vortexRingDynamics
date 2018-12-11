"""
This is the main script for ASNE 5051 Project 2. It implements a simple
Biot-Savart Law solver to analyze the interaction between two vortex rings.

Currently, this solver exploits the symmetry that is known to exsist in the case
of perfectly circular rings - although this is stargithforward to generalize further
in the future.

Authors: Lucas Calvert, Duncan McGough

Created: 12/3/18
Edited: 12/11/18
"""

######################################
# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import trapz
import time
from matplotlib.widgets import Slider, Button, RadioButtons
######################################



######################################
# GIVEN INFO, CONSTANTS AND OTHER VARIABLES


#These are the input parameters that may change
H0 = 1 #initial spacing of the rings
ring1dir = 1 # direction of rotation of vortex ring one
ring2dir = 1 # direction of rotation of vortex ring one #MAKE NEGATIVE to change problem setup
R0 = 1 #initial radius of both rings
Gamma = 1.0 #vortex strength for both rings
tStep = 0.1 #time step for each iteration
displayValue = 1 #The plot will update every "displayValue" iterations
res = 100 #nubmer of points used to model each vortex ring
tTot = 25.0 #total time of simulation


#Other important definitions
pi = np.pi #pi
numTime = np.int(tTot/tStep)
dim = 3 #3D problem
initialSleepValue = 0 # [seconds] Set this to pause plotting and simulation on first iteration to help get a better screen capture

#Pre-allocate arrays for both rings. The 3D part is for each time step. Rows are x,y,z respectively, and each column corresponds to the points of the ring
ring1 = np.zeros((numTime,3,res))
ring2 = np.zeros((numTime,3,res))


#note that that thetaP is a single location, and we integrate about theta in thi code.
theta = np.linspace(0.0,2*pi,res) #defining a vector of angles to make a circle
thetaP = 0.0 #lets just call this zero since it shouldn't make a difference
######################################




######################################
# FUNCTIONS
def biotSavart(R,Rp,H,Hp,ringdir,ringdirp,deltaT,Gamma):


    x = np.array([np.ones(res)*H,R*np.cos(theta),R*np.sin(theta)])

    # Note that, since we can exploit symmetry, we only have to calculate everythin at a singe values of thetaPrime
    xP = np.array([np.ones(res)*Hp,np.ones(res)*Rp*np.cos(thetaP),np.ones(res)*Rp*np.sin(thetaP)])

    #Direction of vorticity vector at x'
    eW = ringdirp*np.array([np.zeros(res),np.ones(res)*np.sin(thetaP),-np.ones(res)*np.cos(thetaP)])
    #eW = ringdirp*np.array([np.zeros(res),np.sin(theta),-np.cos(theta)])


    integrand = np.zeros((3,res))
    #Lets setup the integrand for biot-savart
    for i in range(res):
        if (np.linalg.norm((x[:,i]-xP[:,i]))**3 == 0):
            integrand[:,i] = 0
        else:
            integrand[:,i] = R*np.divide(np.cross(eW[:,i],(x[:,i]-xP[:,i])) , np.linalg.norm((x[:,i]-xP[:,i]))**3 )



    #Now, lets estimate the value of the biot-savart integral
    Vinduced = np.divide(Gamma,(4*pi))*trapz(integrand, dx=(2*pi/res))
    #print(Vinduced)


    # Find change in radius and position of the ring that the velocity is being induced AT
    deltaR = Vinduced[1]*deltaT
    deltaH = Vinduced[0]*deltaT


    # Return the results
    return deltaR, deltaH



######################################






#NOTE:the vortices are both defined in the y-z plane, and "leapfrog" in the x direction

fig1 = plt.figure(figsize=(20,5))
ringView = fig1.add_subplot(133, projection='3d')
radiusView = fig1.add_subplot(132)
distView = fig1.add_subplot(131)

#Lets print the initial conditions on the left side of all of the plots
bottomLeftY = 0.3
bottomLeftX = 0.01
boxHeight = 0.04
boxWidth = 0.07
vSpacing = 0.01
numButtons = 7
startX = np.ones(numButtons)
startY = np.ones(numButtons)
print(startX)
for i in range(numButtons):
    print(i)
    startX[i] = bottomLeftX
    startY[i] = bottomLeftY + (boxHeight + vSpacing)*i

H0BoxAx = plt.axes([startX[0],startY[0], boxWidth , boxHeight])
R0BoxAx = plt.axes([startX[1],startY[1], boxWidth , boxHeight])
Gamma1BoxAx = plt.axes([startX[2],startY[2], boxWidth , boxHeight])
Gamma2BoxAx = plt.axes([startX[3],startY[3], boxWidth , boxHeight])
tStepBoxAx = plt.axes([startX[4],startY[4], boxWidth , boxHeight])
tTotBoxAx = plt.axes([startX[5],startY[5], boxWidth , boxHeight])
resBoxAx = plt.axes([startX[6],startY[6], boxWidth , boxHeight])

H0BoxButton = Button(H0BoxAx, "H0 = " + str(H0))
R0BoxButton = Button(R0BoxAx, "R0 = " + str(R0))
Gamma1Button = Button(Gamma1BoxAx, "Gamma 1 = " + str(Gamma*ring1dir))
Gamma2Button = Button(Gamma2BoxAx, "Gamma 2 = " + str(Gamma*ring2dir))
tStepButton = Button(tStepBoxAx, "t step = " + str(tStep))
tTotButton = Button(tTotBoxAx, "Total t = " + str(tTot))
resBoxButton= Button(resBoxAx, "# Points = " + str(res))


######################################
## Set axis INFO

#ringView
ringView.set_title('Vortex Rings in Space')

ringView.set_xlim3d([-1, 20])
ringView.set_xlabel('X [m]')

ringView.set_ylim3d([-2, 2])
ringView.set_ylabel('Y [m]')

ringView.set_zlim3d([-2, 2])
ringView.set_zlabel('Z [m]')

#distance between rings view
distView.set_title('Horizontal Distance Between Rings')

distView.set_xlim([0, tTot])
distView.set_xlabel('time [s]')

distView.set_ylim([-3*H0, 3*H0])
distView.set_ylabel('Distance [m]')


#radius of each ring view

radiusView.set_title('Radius of Each Ring')

radiusView.set_xlim([0, tTot])
radiusView.set_xlabel('time [s]')

radiusView.set_ylim([0, 2*R0])
radiusView.set_ylabel('Radius [m]')


######################################

R1 = R0
R2 = R0
H1 = 0.0
H2 = H0
tCurrent = 0.0



for i in range(numTime):

    deltaT = tStep
    tCurrent = tCurrent+deltaT

    #Calculate changes in height and radius due to all 4 induced velocities
    deltaR1,deltaH1 = biotSavart(R1,R1,H1,H1,ring1dir,ring1dir,deltaT,Gamma)
    deltaR2,deltaH2 = biotSavart(R1,R2,H1,H2,ring1dir,ring2dir,deltaT,Gamma)
    deltaR3,deltaH3 = biotSavart(R2,R1,H2,H1,ring2dir,ring1dir,deltaT,Gamma)
    deltaR4,deltaH4 = biotSavart(R2,R2,H2,H2,ring2dir,ring2dir,deltaT,Gamma)

    #Calculating new parameters of the first ring
    R1 = R1 - (deltaR1+deltaR2)
    H1 = H1 - (deltaH1+deltaH2)
    #Calculating new parameter of the second ring
    R2 = R2 - (deltaR3+deltaR4)
    H2 = H2 - (deltaH3+deltaH4)

    currentDist = H2-H1

    #Update the plot every "displayValue" calculation iterations
    if (np.mod(i,displayValue) == 0):



        ## Find coords of ring1
        y1 =  R1*np.cos(theta) #x values of points for the first ring (BACK)
        z1 =  R1*np.sin(theta) #y values of points for the first ring (BACK)
        x1 = np.zeros(res) + H1

        #Find coords of ring2
        y2 =  R2*np.cos(theta) #x values of points for the second ring (FRONT)
        z2 =  R2*np.sin(theta) #y values of points for the second ring (FRONT)
        x2 = np.zeros(res) + H2

        ## Display updated time
        ctrString = "Time: " + str(tCurrent) + "s"
        counterText = radiusView.text(0.05, 0.93,ctrString, transform=radiusView.transAxes)


        ## PLOTTING
        scat1 = ringView.scatter(x1,y1,z1,c="b")
        scat2 = ringView.scatter(x2,y2,z2,c="r")
        distScat = distView.scatter(tCurrent,currentDist,c='k')
        radScat1 = radiusView.scatter(tCurrent,R1,c="b")
        radScat2 = radiusView.scatter(tCurrent,R2,c="r")

        plt.pause(0.001)

        #if (R1<=1.01) and (R1>=0.99):
        #    if (R2<=1.01 and R2>=0.99):
        #        print("Leapfrog Cycle Complete at t =",tCurrent)

        if ring2dir == 1:
            if (np.abs(R1-R2)<=0.05):
                print("Leapfrog Cycle Complete at t =",tCurrent)
        else:
            if (np.abs(currentDist - H0/2)<=0.02):
                print("H = H/2 at t =",tCurrent)


        if (i==0):
            time.sleep(initialSleepValue)

        scat1.remove()
        scat2.remove()
        counterText.remove()


#Hold the display ON after simulation ends
plt.show()
