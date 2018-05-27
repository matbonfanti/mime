
import sys

# Scipy modules
from sympy import *
import numpy as np
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

# interface to fitting functions
from optimization import *

#######################################################################################################
class sumofproduct:
   """ define and fit 2D set of data with a sum of product analytical form """
#######################################################################################################

   def __init__(self, xBasisLabel, yBasisLabel ):

      # Define symbols for the variable
      self.x = symbols("x")
      self.y = symbols("y")

      self.xBasis = []
      self.yBasis = []
      for label in xBasisLabel:
         self.xBasis.append( sympify(label) )
      for label in yBasisLabel:
         self.yBasis.append( sympify(label) )

      self.fullBasis = []
      for fx in self.xBasis:
         for fy in self.yBasis:
            self.fullBasis.append( fx*fy )


   def fit_data(self, Coord, V ):

      # check that the number of points is equal to the number of potential value
      if len(Coord) != len(V):
         print(" Wrong number of data points... "); sys.exit()
      # check that all the points have two dimensional coordinate x and y
      for point in Coord:
         if len(point)  != 2:
            print(" Wrong coordinate dimension... "); sys.exit()

      # do the fitting
      return linearfitting( Coord, V, self.fullBasis, (self.x, self.y) )


   def symbolic_expr(self, Params):

      # define global potential function
      expr = 0.
      for val, fun in zip(Params,self.fullBasis):
         expr += val * fun

      # return the symbolic expression of the fitting function
      return expr


   def plot_function(self, Params, xLim = (0.,1.), yLim = (0.,1.), nx = 10, ny = 10 ):

      # create symbolic expression of the fitting function
      expr = self.symbolic_expr( Params )
      # create numpy function from symbolic expression
      func = lambdify( (self.x, self.y), expr, "numpy" )

      # define data grid
      xgrid = np.linspace( xLim[0], xLim[1] , nx )
      ygrid = np.linspace( yLim[0], yLim[1] , ny )
      X, Y = np.meshgrid(xgrid, ygrid)
      V = func(X, Y)

      # create plot
      fig, ax = plt.subplots(1, 1)
      CS = ax.contour(X, Y, V, 10, colors='black')
      plt.clabel(CS, inline=1, fontsize=10)
      plt.contourf(X, Y, V, 100, cmap='RdBu_r')
      plt.colorbar();

      # save figure to file
      plt.savefig( "plot_fitting_function.pdf", format="pdf" )

      #plt.show()


   def plot_residue(self, Coord, Pot, Params, xLim = (0.,1.), yLim = (0.,1.), nx = 10, ny = 10 ):

      # create symbolic expression of the fitting function
      expr = self.symbolic_expr( Params )
      # create numpy function from symbolic expression
      func = lambdify( (self.x, self.y), expr, "numpy" )

      # define data grid
      X = []; Y = []; V = []
      for x,v in zip(Coord,Pot):
         X.append(x[0])
         Y.append(x[1])
         V.append( func(x[0], x[1]) - v )

      spline_res = interp2d(X, Y, V)

      # define data grid
      xgrid = np.linspace( xLim[0], xLim[1] , nx )
      ygrid = np.linspace( yLim[0], yLim[1] , ny )
      V = []
      for y in ygrid:
         V.append( [ float(spline_res(x, y)) for x in xgrid ] )
      V = np.array(V)

      # create plot
      fig, ax = plt.subplots(1, 1)
      CS = ax.contour(xgrid, ygrid, V, 10, colors='black')
      plt.clabel(CS, inline=1, fontsize=10)
      plt.contourf(xgrid, ygrid, V, 100, cmap='RdBu_r')
      plt.colorbar();

      # save figure to file
      plt.savefig( "plot_residue.pdf", format="pdf" )
