
from sympy import *
import numpy as np
from numpy.linalg import lstsq

#*****************************************************************************

def linearfitting( X, V, BasisFuncExpr, Var, NonLinearValues = {}, Weights = [] ):
   """ Computes the least-square optimal linear combinations
      of the functions BasisFuncExpr (defined with sympy symbolic expr)
      of the variable Var (also defined with its symbol),
      fitting the points ( X(i), V(i) ) defined in the lists X and V.
      The dictionary NonLinearValues defines the values of the eventual
      non linear parameters as { symbol:value } where symbol is the
      symbolic expression of the parameters and value its fixed value """

   # Set the vector of the weights
   if len(Weights) < len(X):
       Weights = [ 1.0 ]*len(X)
   else:
       Weights = np.array(Weights)*len(X)/sum(np.array(Weights))

   # Set the numpy functions of the basis set functions
   BasisFunc = []
   for Basis in BasisFuncExpr:

      # substitute the values of the non linear parameters
      LinExpr = Basis
      for NonLinPar in NonLinearValues.keys():
         LinExpr = LinExpr.subs( NonLinPar, NonLinearValues[NonLinPar] )
      # Transform the expression in a numpy function
      BasisFunc.append( lambdify( Var, LinExpr, "numpy" ) )

   # Define the matrix of the coefficients
   A = []
   for f in BasisFunc:
      if isinstance( Var, list ) or isinstance( Var, tuple ):
         A.append( [ w*f( *singleX ) for singleX,w in zip(X,Weights) ] )
      else:
         A.append( [ f( singleX ) for singleX,w in zip(X,Weights) ] )
   A = np.array( A ).transpose()

   # Define the vector of the V values
   B = np.array( [v*w for v,w in zip(V,Weights)] )

   # Solve the AX = B system with least-squares and return the solution
   Coefficients, ChiSquared, Rank, Singular = lstsq( A, B )

   return Coefficients, ChiSquared[0]


#*****************************************************************************

def nonlinearfitting( X, V, FittingFunction, Var, P0, ConstrainedParameters = {}, Weights = [] ):
   """ Computes the least-square optimal linear combinations
      of the functions BasisFuncExpr (defined with sympy symbolic expr)
      of the variable Var (also defined with its symbol),
      fitting the points ( X(i), V(i) ) defined in the lists X and V.
      The dictionary NonLinearValues defines the values of the eventual
      non linear parameters as { symbol:value } where symbol is the
      symbolic expression of the parameters and value its fixed value """

   # Set the vector of the weights
   if len(Weights) < len(X):
       Weights = [ 1.0 ]*len(X)
   else:
       Weights = np.array([ 1.0/w**2 for w in Weights ])/len(X)*sum(np.array(Weights)**2)

   # Substitute linear parameters
   for Param in ConstrainedParameters.keys():
         FittingFunction = FittingFunction.subs( Param, ConstrainedParameters[Param] )

   # Transform the expression in a numpy function
   Func = lambdify( Var, FittingFunction, "numpy" )

   # do the fitting with scipy curve_fit
   Coefficients, Covariance = curve_fit( Func, X, V, p0=P0, sigma=Weights )

   return Coefficients, sqrt(diag(Covariance))

#*****************************************************************************

