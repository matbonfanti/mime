#! /usr/bin/env python

import math
import matplotlib.pyplot as plt
import numpy as np

# name of the data file that is created by the function of this module
DATAFILE = "test_data.dat"
DATAPLOT = "test_data.pdf"

# interval of definition of the input data
XLIMITS = (-2., +2.)
YLIMITS = (0., 90.)

# some parameters
_A0 = 0.2
_XCENTER0 = 0.3
_DELTAA = 0.1
_DELTAX = 1.5

# fix random seed
np.random.seed(42)


# =====================================================================================================================

# define an arbitrary function depending on two variables: a distance and an angular value (function has random noise)
# this function will produce the 2D data that will be then fit with mime
def potentialfunction(x0, theta):
    a = _A0 + _DELTAA * math.cos(theta / 180. * math.pi)
    xcenter = _XCENTER0 + _DELTAX * math.sin(theta / 180. * math.pi)
    v = a * (x0 - xcenter) ** 2 + np.random.random_sample() * 0.05
    return v


# =====================================================================================================================

# define the function that will produce the data,
def writedatafile():

    # define a uniform 2D grid
    xgrid = np.linspace(XLIMITS[0], XLIMITS[1], 20)
    ygrid = np.linspace(YLIMITS[0], YLIMITS[1], 18)
    X, Y = np.meshgrid(xgrid, ygrid)

    # vectorize the function with numpy
    potentialvec = np.vectorize(potentialfunction)
    # evaluate the function on the grid
    V = potentialvec(X, Y)

    # write the data to file, for later use
    with open(DATAFILE, "w") as f:
        for xln, yln, vln in zip(X, Y, V):
            for x, y, z in zip(xln, yln, vln):
                f.write("%15.7f%15.7f%15.7f\n" % (x, y, z))

    # create a 2d plot of the function and show the plot to screen
    fig, ax = plt.subplots(1, 1)
    CS = ax.contour(X, Y, V, 10, colors='black')
    plt.clabel(CS, inline=1, fontsize=10)
    plt.contourf(X, Y, V, 100, cmap='RdBu_r')
    plt.colorbar()
    plt.savefig(DATAPLOT, format="pdf")


# =====================================================================================================================

if __name__ == '__main__':
    writedatafile()
