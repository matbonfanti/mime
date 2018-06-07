
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


def write_mctdh_opfile( nMonomer, nBath, GSfunc, ESfunc, COUPfunc, RingGS_Morse, RingES_Morse ):

   # prepare 1D and 2D functions needed to construct the hamiltonian

   torsion_eq = 0.0

   # 2D function for ground state potential as a function of torsion and stretching
   GS_2D = GSfunc.symbolic_expr( ).as_coefficients_dict()
   # 1D function for the ground state potential
   GS_1D_x = GSfunc.symbolic_expr( ).subs(GSfunc.y,torsion_eq).as_coefficients_dict()

   # 2D function for excited state potential as a function of torsion and stretching
   ES_2D = ESfunc.symbolic_expr( ).as_coefficients_dict()
   # 1D function for the excited state potential as a function of stretching
   ES_1D_x = ESfunc.symbolic_expr( ).subs(ESfunc.y,torsion_eq).as_coefficients_dict()

   # 2D function for coupling potential as a function of torsion and stretching
   COUP_2D = COUPfunc.symbolic_expr( ).as_coefficients_dict()
   # 1D function for the coupling potential as a function of stretching
   COUP_1D_x = COUPfunc.symbolic_expr( ).subs(COUPfunc.y,torsion_eq).as_coefficients_dict()

   # define dictionary for potential printing, GROUND STATE
   GS_1D_x_print = []
   for monom,n in zip(GS_1D_x,range(len(GS_1D_x))):
      GS_1D_x_print.append( [ "gs_str_a{}".format(n), GS_1D_x[monom], GSfunc.xBasisFmtDict[monom] ] )
   GS_2D_print = []
   for term,n in zip(GS_2D,range(len(GS_2D))):
      fun_tor, fun_str = term.as_independent( GSfunc.x )
      GS_2D_print.append( [ "gs_strtor_b{}".format(n), GS_2D[term], GSfunc.xBasisFmtDict[fun_str], GSfunc.yBasisFmtDict[fun_tor] ] )

   # define dictionary for potential printing, EXCITED STATE
   ES_1D_x_print = []
   for monom,n in zip(ES_1D_x,range(len(ES_1D_x))):
      ES_1D_x_print.append( [ "es_str_a{}".format(n), ES_1D_x[monom], ESfunc.xBasisFmtDict[monom] ] )
   ES_2D_print = []
   for term,n in zip(ES_2D,range(len(ES_2D))):
      fun_tor, fun_str = term.as_independent( ESfunc.x )
      ES_2D_print.append( [ "es_strtor_b{}".format(n), ES_2D[term], ESfunc.xBasisFmtDict[fun_str], ESfunc.yBasisFmtDict[fun_tor] ] )

   # define dictionary for potential printing, COUPLING
   COUP_1D_x_print = []
   for monom,n in zip(COUP_1D_x,range(len(COUP_1D_x))):
      COUP_1D_x_print.append( [ "coup_str_a{}".format(n), COUP_1D_x[monom], COUPfunc.xBasisFmtDict[monom] ] )
   COUP_2D_print = []
   for term,n in zip(COUP_2D,range(len(COUP_2D))):
      fun_tor, fun_str = term.as_independent( COUPfunc.x )
      COUP_2D_print.append( [ "coup_strtor_b{}".format(n), COUP_2D[term], COUPfunc.xBasisFmtDict[fun_str], COUPfunc.yBasisFmtDict[fun_tor] ] )

   # define function to print lines for 1D
   def stringfor1Dpotential( factor, nr_dof, potential_print, el=0 ):
      string = ""
      for term in potential_print:
         if el == 0:
            string += "{}*{} |{} {} \n".format(factor, term[0], nr_dof, term[2])
         else:
            string += "{}*{} |1 S{}&{} |{} {} \n".format(factor, term[0], el[0], el[1], nr_dof, term[2])
      return string
   # define function to print lines for 2D
   def stringfor2Dpotential( factor, nr_dof1, nr_dof2, potential_print, el=0 ):
      string = ""
      for term in potential_print:
         if el == 0:
            string += "{}*{} |{} {} |{} {} \n".format(factor, term[0], nr_dof1, term[2], nr_dof2, term[3])
         else:
            string += "{}*{} |1 S{}&{} |{} {} |{} {} \n".format(factor, term[0], el[0], el[1], nr_dof1, term[2], nr_dof2, term[3])
      return string


   f = open( str(nMonomer)+"-mer.op", "w" )

   # list of the torsions that are not fixed in the dynamics (nr of the left monomer)
   MovingTorsion = [ nMonomer/2 ]

   # define names of the variables
   elcoord      = "el"
   act_torsions = [ "tor"+"%02d%02d"%(i,i+1) for i in MovingTorsion ]
   b_stretch    = [ "bs"+"%02d%02d"%(i+1,i+2) for i in range(nMonomer-1) ]
   internal     = [ "int"+ "%02d"%(i+1) for i in range(nMonomer) ]
   bath_modes   = [ "bath"+ "%02d"%(i+1) for i in range(nBath*len(MovingTorsion)) ]

   # define numbers of the variables
   nr_elcoord   = 1
   nr_act_torsions = [ nr_elcoord+1+i for i in range(len(MovingTorsion)) ]
   nr_b_stretch    = [ nr_act_torsions[-1] + 1 + i for i in range(len(b_stretch)) ]
   nr_internal     = [ nr_b_stretch[-1] + 1 + i for i in range(len(internal)) ]
   nr_bath_modes   = [ nr_internal[-1] + 1 + i for i in range(len(bath_modes)) ]

   ## ground and excited states reference energy
   #shifttorsgs = -0.0001649497
   #shifttorses = 0.0051169274
   #coupling = -0.0353873606

   f.write( "OP_DEFINE-section\ntitle\n1stack of {}mers\nend-title\nend-op_define-section \n""".format(nMonomer) )

   f.write( "\nPARAMETER-SECTION\n" )

   #f.write( "\n# potential shifts\n" )
   #f.write( "shifttorsgs = {}\n".format(shifttorsgs) )
   #f.write( "shifttorses = {}\n".format(shifttorses) )
   #f.write( "w = {}\n".format(coupling) )

   f.write( "\n# masses of the torsions\n" )
   for label in act_torsions:
      f.write( "mass_{} = 1365700.43360989".format(label) + "\n" )
   f.write( "\n# masses of the bond stretching modes \n" )
   for label in b_stretch:
      f.write( "mass_{} = 13985.5445411623".format(label) + "\n" )
   f.write( "\n# masses of the internal modes \n" )
   for label in internal:
      f.write( "mass_{} = 94846.9644229036".format(label) + "\n" )
   f.write( "\n# constant for the kinetic coupling between stretching and internal modes \n" )
   f.write( "gbsint = -0.0000000109\n" )

   f.write( "\n# coefficients for the ground state potential as a function of the stretching mode \n" )
   for term in GS_1D_x_print:
      f.write( "{} = {}\n".format(term[0], term[1]) )
   f.write( "\n# coefficients for the ground state potential as a function of the stretching mode and torsion mode \n" )
   for term in GS_2D_print:
      f.write( "{} = {}\n".format(term[0], term[1]) )

   f.write( "\n# coefficients for the excited state potential as a function of the stretching mode \n" )
   for term in ES_1D_x_print:
      f.write( "{} = {}\n".format(term[0], term[1]) )
   f.write( "\n# coefficients for the excited state potential as a function of the stretching mode and torsion mode \n" )
   for term in ES_2D_print:
      f.write( "{} = {}\n".format(term[0], term[1]) )

   f.write( "\n# coefficients for the coupling potential as a function of the stretching mode \n" )
   for term in COUP_1D_x_print:
      f.write( "{} = {}\n".format(term[0], term[1]) )
   f.write( "\n# coefficients for the coupling potential as a function of the stretching mode and torsion mode \n" )
   for term in COUP_2D_print:
      f.write( "{} = {}\n".format(term[0], term[1]) )
      
   f.write( "\n# coefficients for the ground state potential as a function of the ring breathing mode \n" )
   f.write( "{} = {}\n".format("digs",     RingGS_Morse[0]) )
   f.write( "{} = {}\n".format("alphaigs", RingGS_Morse[1]) )
   f.write( "{} = {}\n".format("x0igs",    RingGS_Morse[2]) )
   f.write( "{} = {}\n".format("e0igs",    RingGS_Morse[3]) )
   f.write( "\n# coefficients for the excited state potential as a function of the ring breathing mode \n" )
   f.write( "{} = {}\n".format("dies",     RingES_Morse[0]) )
   f.write( "{} = {}\n".format("alphaies", RingES_Morse[1]) )
   f.write( "{} = {}\n".format("x0ies",    RingES_Morse[2]) )
   f.write( "{} = {}\n".format("e0ies",    RingES_Morse[3]) )      
      
   f.write( "\n# parameters for the definition of the bath modes \n" )
   f.write( "mbath = 1.0\n" )
   f.write( "deltaom = 0.00003\n" )
   f.write( "cfac = 0.0355222904\n" )
   f.write( "shi = cfac*cfac/mbath\n" )

   f.write( "\n# frequencies of the bath modes \n" )
   for i in range(len(MovingTorsion)):
      for j in range(nBath):
         f.write( "om{:02d} = {}*deltaom\n".format(i*nBath+j+1, float(j+1)) )
   f.write( "\n# masses of the bath modes \n" )
   for label in bath_modes:
      f.write( "mass_{} = mbath".format(label) + "\n" )
   f.write( "\n# quadratic potential terms of the bath modes \n" )
   for i in range(len(MovingTorsion)*nBath):
      f.write( "k{:02d} = 0.5*mbath*om{:02d}*om{:02d}\n".format(i+1, i+1, i+1) )
   f.write( "\n# coupling coefficients to the system \n" )
   for i in range(len(MovingTorsion)*nBath):
      f.write( "c{:02d} = cfac*om{:02d}\n".format(i+1, i+1) )

   f.write( "\nend-parameter-section\n" )

   f.write( """\nLABELS-SECTION

# morse potentials for the ground and excited state momomer V as a function of the ring breathing internal mode
intgs = morse1[digs,alphaigs,x0igs,e0igs]
intes = morse1[dies,alphaies,x0ies,e0ies]

end-labels-section
""" )

   f.write( "\nHAMILTONIAN-SECTION\n" )
   f.write( "\n# electronic mode\n" )
   f.write( "modes | {}".format(elcoord) + "\n" )

   def write_formatted_names_of_modes( listofnames, size ):
      for sublist in [ listofnames[i:i+size] for i in range(0, len(listofnames), size) ]:
         f.write( "modes "+" ".join( ["| {:^7}".format(s) for s in sublist ] )+"\n" )

   f.write( "# name for the central torsion mode\n" )
   write_formatted_names_of_modes( act_torsions, 10 )
   f.write( "# names for the bond stretching modes\n" )
   write_formatted_names_of_modes( b_stretch, 10 )
   f.write( "# names for the internal monomer modes\n" )
   write_formatted_names_of_modes( internal, 10 )
   f.write( "# names for the bath modes coupled to central torsion\n" )
   write_formatted_names_of_modes( bath_modes, 10 )

   f.write( "\n# kinetic coupling between stretching and internal modes\n" )
   for i in range(nMonomer-1):
      f.write( "-gbsint |{} dq |{} dq\n".format(nr_b_stretch[i],nr_internal[i]) )
      f.write( "-gbsint |{} dq |{} dq\n".format(nr_b_stretch[i],nr_internal[i+1]) )

   f.write( "\n# Kinetic operator for the torsion\n" )
   for i in nr_act_torsions:
      f.write( "1.0 |{} KE\n".format(i) )
   f.write( "\n# Kinetic operator for the stretching coordinates\n" )
   for i in nr_b_stretch:
      f.write( "1.0 |{} KE\n".format(i) )
   f.write( "\n# Kinetic operator for the internal monomer coordinates\n" )
   for i in nr_internal:
      f.write( "1.0 |{} KE\n".format(i) )

   f.write( "\n# Ground state potential\n" )
   #f.write( "# potential shift\n" )
   #f.write( "{}*shifttorsgs |1 1\n".format( 38.0 ) )
   f.write( "# potential of the internal modes\n" )
   for i in range(nMonomer):
      f.write( "1.0 |{} {}\n".format(nr_internal[i], "intgs") )
   f.write( "# potential of the inter-monomer stretching and active torsion modes\n" )
   for i in range(nMonomer-1):
      if i+1 in MovingTorsion:
         tor_i = MovingTorsion.index(i+1)
         f.write( stringfor2Dpotential( 2.0, nr_b_stretch[i], nr_act_torsions[tor_i], GS_2D_print ) )
      else:
         f.write( stringfor1Dpotential( 2.0, nr_b_stretch[i], GS_1D_x_print ) )

   for i in range(nMonomer):
      f.write( "\n# Diagonal excited state potential, exciton on the monomer {}\n".format(i+1) )
      #f.write( "# potential shift\n")
      #f.write( " 2.0*shifttorses |1 S{}&{}\n".format(i+1,i+1) )
      #f.write( "-2.0*shifttorsgs |1 S{}&{}\n".format(i+1,i+1) )
      f.write( "# potential of the internal modes\n" )
      f.write( " 1.0 |1 S{}&{} |{} {} \n".format(i+1,i+1,nr_internal[i], "intes") )
      f.write( "-1.0 |1 S{}&{} |{} {} \n".format(i+1,i+1,nr_internal[i], "intgs") )
      if i != 0:                     # this is not the first monomer, left torsion and stretching do exist
         f.write( "# potential of the left intermonomer modes\n" )
         if i in MovingTorsion:
            tor_i = MovingTorsion.index(i)
            f.write( stringfor2Dpotential(  1.0, nr_b_stretch[i-1], nr_act_torsions[tor_i], ES_2D_print, [i+1,i+1] ) )
            f.write( stringfor2Dpotential( -1.0, nr_b_stretch[i-1], nr_act_torsions[tor_i], GS_2D_print, [i+1,i+1] ) )
         else:
            f.write( stringfor1Dpotential(  1.0, nr_b_stretch[i-1], ES_1D_x_print, [i+1,i+1] ) )
            f.write( stringfor1Dpotential( -1.0, nr_b_stretch[i-1], GS_1D_x_print, [i+1,i+1] ) )
      if i != nMonomer-1:            # this is not the last monomer, right torsion and stretching do exist
         f.write( "# potential of the right intermonomer modes\n" )
         if i+1 in MovingTorsion:
            tor_i = MovingTorsion.index(i+1)
            f.write( stringfor2Dpotential(  1.0, nr_b_stretch[i], nr_act_torsions[tor_i], ES_2D_print, [i+1,i+1] ) )
            f.write( stringfor2Dpotential( -1.0, nr_b_stretch[i], nr_act_torsions[tor_i], GS_2D_print, [i+1,i+1] ) )
         else:
            f.write( stringfor1Dpotential(  1.0, nr_b_stretch[i], ES_1D_x_print, [i+1,i+1] ) )
            f.write( stringfor1Dpotential( -1.0, nr_b_stretch[i], GS_1D_x_print, [i+1,i+1] ) )

   for i in range(nMonomer-1):
      f.write( "\n# Coupling potential between excitons on monomer {} and {}\n".format(i+1,i+2) )
      #f.write( "# coupling\n")
      #f.write( " w |1 S{}&{}\n".format(i+1,i+2) )
      f.write( "# potential of the intermonomer modes\n" )
      if i+1 in MovingTorsion:
         tor_i = MovingTorsion.index(i+1)
         f.write( stringfor2Dpotential( 1.0, nr_b_stretch[i], nr_act_torsions[tor_i], COUP_2D_print, [i+1,i+2] ) )
      else:
         f.write( stringfor1Dpotential( 1.0, nr_b_stretch[i], COUP_1D_x_print, [i+1,i+2] ) )

   f.write( "\n# Kinetic energy of the bath coupled to the torsion of mode {}\n".format(1) )
   for i in nr_bath_modes:
      f.write( "1.0 |{} KE\n".format(i) )

   f.write( "\n# Potential of the bath coupled to the torsion of mode {}\n".format(1) )
   for i,nr in zip(range(len(MovingTorsion)*nBath),nr_bath_modes):
      f.write( "k{:02d} |{} q^2\n".format(nr,i) )
   for i in range(len(MovingTorsion)):
      for j in range(nBath):
         f.write( "-c{:02d} |{} q |{} q \n".format(i*nBath+j+1, nr_act_torsions[i], nr_bath_modes[i*nBath+j]) )
   for i in range(len(MovingTorsion)):
      f.write("20.0*0.5*shi   |{} q^2 \n".format(nr_act_torsions[i]))


   f.write( "\nend-hamiltonian-section\n" )
   f.write( "\nend-operator\n" )

   f.close()
