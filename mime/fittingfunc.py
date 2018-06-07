
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

   def __init__(self, xBasisLabel, yBasisLabel, xBasisFormat = [], yBasisFormat = [] ):

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

      self.xBasisFmtDict = {}
      if len(xBasisFormat) == 0:
         for label,fx in zip(xBasisLabel,self.xBasis):
            self.xBasisFmtDict[fx] = label
      else:
         for label,fx in zip(xBasisFormat,self.xBasis):
            self.xBasisFmtDict[fx] = label

      self.yBasisFmtDict = {}
      if len(yBasisFormat) == 0:
         for label,fy in zip(yBasisLabel,self.yBasis):
            self.yBasisFmtDict[fy] = label
      else:
         for label,fy in zip(yBasisFormat,self.yBasis):
            self.yBasisFmtDict[fy] = label


   def fit_data(self, Coord, V ):

      # check that the number of points is equal to the number of potential value
      if len(Coord) != len(V):
         print(" Wrong number of data points... "); sys.exit()
      # check that all the points have two dimensional coordinate x and y
      for point in Coord:
         if len(point)  != 2:
            print(" Wrong coordinate dimension... "); sys.exit()

      # do the fitting
      self.Coeffs, Chisquared = linearfitting( Coord, V, self.fullBasis, (self.x, self.y) )

      return Chisquared


   def symbolic_expr(self):

      # define global potential function
      expr = 0.
      for val, fun in zip(self.Coeffs,self.fullBasis):
         expr += val * fun

      # return the symbolic expression of the fitting function
      return expr


   def show_1D_modelplusdata(self, Coord, Pot, yspacing=1.0 ):

      # define data grid
      Qdata = []; Vdata = []
      for x,v in zip(Coord,Pot):
         Qdata.append(x[0]+x[1]*yspacing)
         Vdata.append( v ) 

      # create symbolic expression of the fitting function
      expr = self.symbolic_expr( )
      # create numpy function from symbolic expression
      func = lambdify( (self.x, self.y), expr, "numpy" )

      # define model grid
      Qmodel = []; Vmodel = []
      for x,y in Coord:
         Qmodel.append( x+y*yspacing )
         Vmodel.append( func(x, y) ) 
            
      # create plot
      fig, graph = plt.subplots(1, 1)
      graph.plot( Qdata, Vdata, label="input data", ls="", marker="o"  )
      graph.plot( Qmodel, Vmodel, label="fitted model" )
      
      # save figure to file
      plt.show()
      plt.close()


   def plot_function(self, xLim = (0.,1.), yLim = (0.,1.), nx = 10, ny = 10, pdfname = "plot_fitting_function.pdf", 
                           scale_x=1.0, scale_y=1.0, scale_v=1.0 ):

      # create symbolic expression of the fitting function
      expr = self.symbolic_expr( )
      # create numpy function from symbolic expression
      func = lambdify( (self.x, self.y), expr, "numpy" )

      # define data grid
      xgrid = np.linspace( xLim[0], xLim[1] , nx )
      ygrid = np.linspace( yLim[0], yLim[1] , ny )
      X, Y = np.meshgrid(xgrid, ygrid)
      V = func(X, Y)

      # create plot
      fig, ax = plt.subplots(1, 1)
      CS = ax.contour(X*scale_x, Y*scale_y, V*scale_v, 10, colors='black')
      plt.clabel(CS, inline=1, fontsize=10)
      plt.contourf(X*scale_x, Y*scale_y, V*scale_v, 100, cmap='RdBu_r')
      plt.colorbar();

      # save figure to file
      plt.savefig( pdfname, format="pdf" )
      plt.close()


   def plot_residue(self, Coord, Pot, nx = 10, ny = 10, pdfname = "plot_residue.pdf", scale_x=1.0, scale_y=1.0, scale_v=1.0 ):

      # create symbolic expression of the fitting function
      expr = self.symbolic_expr( )
      # create numpy function from symbolic expression
      func = lambdify( (self.x, self.y), expr, "numpy" )

      # define data grid
      X = []; Y = []; V = []
      for x,v in zip(Coord,Pot):
         X.append(x[0])
         Y.append(x[1])
         V.append( v- func(x[0], x[1]) )

      spline_res = interp2d(X, Y, V, kind='linear')

      # define data grid
      xgrid = np.linspace( np.amin(X), np.amax(X), nx )
      ygrid = np.linspace( np.amin(Y), np.amax(Y), ny )
      V = []
      for y in ygrid:
         V.append( [ float(spline_res(x, y)) for x in xgrid ] )
      V = np.array(V)

      # create plot
      fig, ax = plt.subplots(1, 1)
      CS = ax.contour(xgrid*scale_x, ygrid*scale_y, V*scale_v, 10, colors='black')
      plt.clabel(CS, inline=1, fontsize=10)
      plt.contourf(xgrid*scale_x, ygrid*scale_y, V*scale_v, 100, cmap='RdBu_r')
      plt.colorbar();

      # save figure to file
      plt.savefig( pdfname, format="pdf" )
      plt.close()

   def get_basis_funcs( self ):
      return self.fullBasis

   def get_coeffs( self ):
      return self.Coeffs


