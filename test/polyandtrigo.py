#! /usr/bin/env python

from math import *
import numpy as np
import mime

# #################################################
#                Input section                    #
# #################################################

# functions for the first degree of freedom (write functions of x)
xFunctionsLabels = ["1", "x**1", "x**2", "x**3"]

# functions for the second degree of freedom (write functions of y)
yFunctionsLabels = ["1", "sin(pi*y/180.)", "cos(pi*y/180.)", "sin(2*pi*y/180.)", "cos(2*pi*y/180.)"]

# file to read the input data
inputfile = "test_data.dat"

# #################################################o

# read and store data from input file
Coord = np.loadtxt(inputfile, usecols=(0, 1))
Pot = np.loadtxt(inputfile, usecols=2)

# define fitting function
fittingf = mime.sumofproduct(xFunctionsLabels, yFunctionsLabels)

# do fitting
Chisquared = fittingf.fit_data(Coord, Pot)
print(" Fitting done... RMS deviation is %f" % sqrt(Chisquared))

# plot fitting function and residue
fittingf.plot_function(xLim=(-2.0, +2.0), yLim=(0.0, 90.0), nx=20, ny=18)
fittingf.plot_residue(Coord, Pot, nx=20, ny=18)
