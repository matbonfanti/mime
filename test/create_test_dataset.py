#! /usr/bin/env python

from math import *
import matplotlib.pyplot as plt
import numpy as np

# parameters
a0       = 0.2
xcenter0 = 0.3
deltaa   = 0.1
deltax   = 1.5

def potentialfunction( x, theta ):
   a       = a0 + deltaa * cos(theta/180.*pi)
   xcenter = xcenter0 + deltax * sin( theta/180.*pi )
   v = a * ( x - xcenter )**2 + np.random.random_sample()*0.05
   return v
potentialvec = np.vectorize( potentialfunction )

xgrid = np.linspace( -2. , +2. , 20 )
ygrid = np.linspace( 0. , 90. , 18 )
X, Y = np.meshgrid(xgrid, ygrid)
V = potentialvec(X, Y)

f = open( "test_data.dat", "w" )
for xln, yln, vln in zip(X,Y,V):
   for x,y,z in zip(xln, yln, vln):
      f.write("%15.7f%15.7f%15.7f\n" % (x,y,z))
f.close()

fig, ax = plt.subplots(1, 1)
CS = ax.contour(X, Y, V, 10, colors='black')
plt.clabel(CS, inline=1, fontsize=10)

plt.contourf(X, Y, V, 100, cmap='RdBu_r')
plt.colorbar();

plt.show()
