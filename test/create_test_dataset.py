#! /usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np

# parameters
a0 = 0.2
xcenter0 = 0.3
deltaa = 0.1
deltax = 1.5


# define an arbitrary function depending on two variables: a distance and an angular value (function has random noise)
def potentialfunction(x0, theta):
    a = a0 + deltaa * math.cos(theta / 180. * math.pi)
    xcenter = xcenter0 + deltax * math.sin(theta / 180. * math.pi)
    v = a * (x0 - xcenter) ** 2 + np.random.random_sample() * 0.05
    return v


# vectorize the function with numpy
potentialvec = np.vectorize(potentialfunction)

# evaluate the function on a uniform 2D grid
xgrid = np.linspace(-2., +2., 20)
ygrid = np.linspace(0., 90., 18)
X, Y = np.meshgrid(xgrid, ygrid)
V = potentialvec(X, Y)

# write the data to file, for later use
f = open("test_data.dat", "w")
for xln, yln, vln in zip(X, Y, V):
    for x, y, z in zip(xln, yln, vln):
        f.write("%15.7f%15.7f%15.7f\n" % (x, y, z))
f.close()

# create a 2d plot of the function and show the plot to screen
fig, ax = plt.subplots(1, 1)
CS = ax.contour(X, Y, V, 10, colors='black')
plt.clabel(CS, inline=1, fontsize=10)
plt.contourf(X, Y, V, 100, cmap='RdBu_r')
plt.colorbar()
plt.show()
