import matplotlib.pyplot as plt
import time
import random



ysample = random.sample(range(-50, 50), 100)

xdata = []
ydata = []

plt.show()

axes = plt.gca()
axes.set_xlim(0, 100)
axes.set_ylim(-50, +50)
line, = axes.plot(xdata, ydata, 'ro')

for i in range(100):
    xdata.append(i)
    ydata.append(ysample[i])
    #xdata = i
    #ydata = ysample[i]
    line.set_xdata(xdata)
    line.set_ydata(ydata)
    plt.draw()
    plt.pause(1e-17)
    time.sleep(0.01)

# add this if you don't want the window to disappear at the end
plt.show()


'''
plt.ion()

for i in range(10):
    scat = plt.scatter(i, i)
    plt.pause(0.5)
    scat.remove()
'''