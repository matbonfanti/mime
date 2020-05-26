#! /usr/bin/env python

import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import mime
from write_mctdh_opfile import *
from sympy import *
from math import pi

interactiveplots = False
pdfplots = True

# #################################################
#              2D V_G potential Fit               #
# #################################################

print("\n *** 2D V_G potential fit ***\n")

# functions for the first degree of freedom (write functions of x)
xFunctLabels = ["1", "x**1", "x**2"]  # , "x**4"  ]
xMCTDHLabels = ["1", "q^1", "q^2"]  # ,  "q^4"  ]

# functions for the second degree of freedom (write functions of y)
yFunctLabels = ["1", "cos(y)", "cos(2*y)", "cos(3*y)", "cos(4*y)"]
yMCTDHLabels = ["1", "cos[1.0,0.0]", "cos[2.0,0.0]", "cos[3.0,0.0]", "cos[4.0,0.0]"]

# define fitting function
VG_2D = mime.sumofproduct(xFunctLabels, yFunctLabels, xBasisFormat=xMCTDHLabels, yBasisFormat=yMCTDHLabels)

# file to read the input data
inputfile = "PPV_intra2D_VG_VE_w.dat"
# read and store data from input file
Coord = np.loadtxt(inputfile, usecols=(1, 0))
Pot = np.loadtxt(inputfile, usecols=2)
# convert angles from deg to rad
Coord = np.array([[x, theta * pi / 180.0] for x, theta in Coord])

# do fitting and get optimal coefficients
Chisquared = VG_2D.fit_data(Coord, Pot, FixZero=False)
NDoF = len(Pot) - len(xFunctLabels) * len(yFunctLabels)
print(" Fitting done... regression standard error is {} eV".format(sqrt(Chisquared / NDoF) * 27.21142))

print("\n{:>15s} | {:60s}".format("coeff.s", "basis functions"))
for c, f in zip(VG_2D.get_coeffs(), VG_2D.get_basis_funcs()):
    print("{:15.6e} | {:60s}".format(c, str(f)))

# visualize data and model 
if interactiveplots:
    VG_2D.show_1D_modelplusdata(Coord, Pot, yspacing=180.0 / pi)

# plot graphs
if pdfplots:
    Name = "VG_2D"
    print("\n Writing potential to file {}".format('"' + Name + '_func.pdf"'))
    VG_2D.plot_function(xLim=(-0.3, +0.3), yLim=(0.0, 0.5 * pi), nx=27, ny=19, pdfname=Name + "_func.pdf",
                        scale_y=180.0 / pi, scale_v=27.21142)
    print(" Writing fitting residues to file {}\n".format('"' + Name + '_residue.pdf"'))
    VG_2D.plot_residue(Coord, Pot, nx=27, ny=19, pdfname=Name + "_residue.pdf", scale_y=180.0 / pi, scale_v=27.21142)

# #################################################o

print("\n *** 2D V_E potential fit ***\n")

# functions for the first degree of freedom (write functions of x)
xFunctLabels = ["1", "x**1", "x**2", "x**3", "x**4"]
xMCTDHLabels = ["1", "q^1", "q^2", "q^3", "q^4"]

# functions for the second degree of freedom (write functions of y)
yFunctLabels = ["1", "cos(y)", "cos(2*y)", "cos(4*y)", "cos(6*y)"]
yMCTDHLabels = ["1", "cos[1.0,0.0]", "cos[2.0,0.0]", "cos[4.0,0.0]", "cos[6.0,0.0]"]

# define fitting function
VE_2D = mime.sumofproduct(xFunctLabels, yFunctLabels, xBasisFormat=xMCTDHLabels, yBasisFormat=yMCTDHLabels)

# file to read the input data
inputfile = "PPV_intra2D_VG_VE_w.dat"
# read and store data from input file
Coord = np.loadtxt(inputfile, usecols=(1, 0))
Pot = np.loadtxt(inputfile, usecols=3)
# convert angles from deg to rad
Coord = np.array([[x, theta * pi / 180.0] for x, theta in Coord])

# do fitting and get optimal coefficients
Chisquared = VE_2D.fit_data(Coord, Pot, FixZero=False)
NDoF = len(Pot) - len(xFunctLabels) * len(yFunctLabels)
print(" Fitting done... regression standard error is {} eV".format(sqrt(Chisquared / NDoF) * 27.21142))

print("\n{:>15s} | {:60s}".format("coeff.s", "basis functions"))
for c, f in zip(VE_2D.get_coeffs(), VE_2D.get_basis_funcs()):
    print("{:15.6e} | {:60s}".format(c, str(f)))

# visualize data and model 
if interactiveplots:
    VE_2D.show_1D_modelplusdata(Coord, Pot, yspacing=180.0 / pi)

# plot graphs
if pdfplots:
    Name = "VE_2D"
    print("\n Writing potential to file {}".format('"' + Name + '_func.pdf"'))
    VE_2D.plot_function(xLim=(-0.3, +0.3), yLim=(0.0, 0.5 * pi), nx=27, ny=19, pdfname=Name + "_func.pdf",
                        scale_y=180.0 / pi, scale_v=27.21142)
    print(" Writing fitting residues to file {}\n".format('"' + Name + '_residue.pdf"'))
    VE_2D.plot_residue(Coord, Pot, nx=27, ny=19, pdfname=Name + "_residue.pdf", scale_y=180.0 / pi, scale_v=27.21142)

# #################################################o

print("\n *** 2D w potential fit ***\n")

# functions for the first degree of freedom (write functions of x)
xFunctLabels = ["1", "x**1", "x**2", "x**3", "x**4"]
xMCTDHLabels = ["1", "q^1", "q^2", "q^3", "q^4"]

# functions for the second degree of freedom (write functions of y)
yFunctLabels = ["1", "cos(y)", "cos(2*y)", "cos(3*y)", "cos(4*y)", "cos(5*y)", "cos(6*y)"]
yMCTDHLabels = ["1", "cos[1.0,0.0]", "cos[2.0,0.0]", "cos[3.0,0.0]", "cos[4.0,0.0]", "cos[5.0,0.0]", "cos[6.0,0.0]"]

# define fitting function
w_2D = mime.sumofproduct(xFunctLabels, yFunctLabels, xBasisFormat=xMCTDHLabels, yBasisFormat=yMCTDHLabels)

# file to read the input data
inputfile = "PPV_intra2D_VG_VE_w.dat"
# read and store data from input file
Coord = np.loadtxt(inputfile, usecols=(1, 0))
Pot = np.loadtxt(inputfile, usecols=4)
# convert angles from deg to rad
Coord = np.array([[x, theta * pi / 180.0] for x, theta in Coord])

# do fitting and get optimal coefficients
Chisquared = w_2D.fit_data(Coord, Pot)
NDoF = len(Pot) - len(xFunctLabels) * len(yFunctLabels)
print(" Fitting done... regression standard error is {} eV".format(sqrt(Chisquared / NDoF) * 27.21142))

print("\n{:>15s} | {:60s}".format("coeff.s", "basis functions"))
for c, f in zip(w_2D.get_coeffs(), w_2D.get_basis_funcs()):
    print("{:15.6e} | {:60s}".format(c, str(f)))

# visualize data and model 
if interactiveplots:
    w_2D.show_1D_modelplusdata(Coord, Pot, yspacing=180.0 / pi)

# plot graphs
if pdfplots:
    Name = "w_2D"
    print("\n Writing potential to file {}".format('"' + Name + '_func.pdf"'))
    w_2D.plot_function(xLim=(-0.3, +0.3), yLim=(0.0, 0.5 * pi), nx=27, ny=19, pdfname=Name + "_func.pdf",
                       scale_y=180.0 / pi, scale_v=27.21142)
    print(" Writing fitting residues to file {}\n".format('"' + Name + '_residue.pdf"'))
    w_2D.plot_residue(Coord, Pot, nx=27, ny=19, pdfname=Name + "_residue.pdf", scale_y=180.0 / pi, scale_v=27.21142)

# #################################################o

print("\n *** 1D ring breathing ground and exc states potential fit ***\n")


def MorsePotential1(Crd, D, alpha):
    return D * (np.exp(-alpha * Crd) - 1.0) ** 2


def MorsePotential2(Crd, D, alpha, x0, E0):
    return D * (np.exp(-alpha * (Crd - x0)) - 1.0) ** 2


# file to read the input data
inputfile = "PPV_ring_VG_VE.dat"
# read and store data from input file
Coord = np.loadtxt(inputfile, usecols=0)
PotGS = np.loadtxt(inputfile, usecols=1)
PotES = np.loadtxt(inputfile, usecols=2)

RingGS_Morse, _ = curve_fit(MorsePotential1, Coord, PotGS, p0=(1.0, 1.0))
RingES_Morse, _ = curve_fit(MorsePotential2, Coord, PotES, p0=(1.0, 1.0, 0.0, 0.0))

if interactiveplots:
    import matplotlib.pyplot as plt

    fig, graph = plt.subplots(1, 1)
    graph.plot(Coord, PotGS, label="ring, GS", ls="", color="r", marker="o")
    graph.plot(Coord, PotES, label="ring, ES", ls="", color="b", marker="o")
    CoordModel = np.linspace(np.amin(Coord), np.amax(Coord), 200)
    ModelGS = MorsePotential1(CoordModel, *RingGS_Morse)
    ModelES = MorsePotential2(CoordModel, *RingES_Morse)
    graph.plot(CoordModel, ModelGS, label="fit GS", color="r")
    graph.plot(CoordModel, ModelES, label="fit ES", color="b")
    plt.show()
    plt.close()

RingGS_Morse = list(RingGS_Morse)
RingGS_Morse.append(0.0), RingGS_Morse.append(0.0)
RingGS_Morse = np.array(RingGS_Morse)

print("\n{:>15s} | {:>15s} | {:>15s}".format("coeff.s", "GS fit", "ES fit"))
print("{:>15s} | {:15.6f} | {:15.6f}".format("D", RingGS_Morse[0], RingES_Morse[0]))
print("{:>15s} | {:15.6f} | {:15.6f}".format("alpha", RingGS_Morse[1], RingES_Morse[1]))
print("{:>15s} | {:15.6f} | {:15.6f}".format("x0", RingGS_Morse[2], RingES_Morse[2]))
print("{:>15s} | {:15.6f} | {:15.6f}".format("E0", RingGS_Morse[3], RingES_Morse[3]))

# #################################################o

print("\n *** MCTDH operator ***\n")

NameOperator = "20-mer.op"
print(" Writing operator file to file {} ".format(NameOperator))

write_mctdh_opfile(20, 20, VG_2D, VE_2D, w_2D, RingGS_Morse, RingES_Morse)

# #################################################o

print("\n *** Frequencies of the 2D torsion-stretching potential ***\n")

fVG = lambdify((VG_2D.x, VG_2D.y), VG_2D.symbolic_expr(), "numpy")
fVE = lambdify((VE_2D.x, VE_2D.y), VE_2D.symbolic_expr(), "numpy")

fVGdx2 = lambdify((VG_2D.x, VG_2D.y), diff(VG_2D.symbolic_expr(), VG_2D.x, VG_2D.x), "numpy")
fVGdy2 = lambdify((VG_2D.x, VG_2D.y), diff(VG_2D.symbolic_expr(), VG_2D.y, VG_2D.y), "numpy")
fVGdxdy = lambdify((VG_2D.x, VG_2D.y), diff(VG_2D.symbolic_expr(), VG_2D.x, VG_2D.y), "numpy")

fVEdx2 = lambdify((VE_2D.x, VE_2D.y), diff(VE_2D.symbolic_expr(), VE_2D.x, VE_2D.x), "numpy")
fVEdy2 = lambdify((VE_2D.x, VE_2D.y), diff(VE_2D.symbolic_expr(), VE_2D.y, VE_2D.y), "numpy")
fVEdxdy = lambdify((VE_2D.x, VE_2D.y), diff(VE_2D.symbolic_expr(), VE_2D.x, VE_2D.y), "numpy")

if interactiveplots:
    # create plot
    import matplotlib.pyplot as plt

    # define data grid
    xgrid = np.linspace(-0.3, +0.3, 27)
    ygrid = np.linspace(-0.5 * pi, 0.5 * pi, 37)
    X, Y = np.meshgrid(xgrid, ygrid)
    V1 = 2 * fVG(X, Y)
    V2 = fVG(X, Y) + fVE(X, Y)

    fig, ax = plt.subplots(1, 1)
    CS = ax.contour(X, Y * 180.0 / pi, V1 * 27.21142, 10, colors='black')
    plt.clabel(CS, inline=1, fontsize=10)
    plt.contourf(X, Y * 180.0 / pi, V1 * 27.21142, 100, cmap='RdBu_r')
    plt.colorbar()
    plt.show()
    plt.close()

    fig2, ax = plt.subplots(1, 1)
    CS = ax.contour(X, Y * 180.0 / pi, V2 * 27.21142, 10, colors='black')
    plt.clabel(CS, inline=1, fontsize=10)
    plt.contourf(X, Y * 180.0 / pi, V2 * 27.21142, 100, cmap='RdBu_r')
    plt.colorbar()
    plt.show()
    plt.close()

mass_stretch = 13985.5445411623
mass_torsion = 1365700.43360989
mixed_mass = sqrt(mass_stretch) * sqrt(mass_torsion)

exc_fraction = np.linspace(0.0, 1.0, 21)
frequency_torsion = []
frequency_stretching = []
for exc_VE_frac in exc_fraction:
    potential = lambda x: (2. - exc_VE_frac) * fVG(*x) + exc_VE_frac * fVE(*x)
    optres = minimize(potential, np.array([0.0, 0.0]))
    print("Minimum of ground state potential is at y = {} bohr and theta = {} deg".format(optres.x[0],
                                                                                          optres.x[1] * 180.0 / pi))
    Hessian = np.array([[((2. - exc_VE_frac) * fVGdx2(*optres.x) + exc_VE_frac * fVEdx2(*optres.x)) / mass_stretch,
                         ((2. - exc_VE_frac) * fVGdxdy(*optres.x) + exc_VE_frac * fVEdxdy(*optres.x)) / mixed_mass],
                        [((2. - exc_VE_frac) * fVGdxdy(*optres.x) + exc_VE_frac * fVEdxdy(*optres.x)) / mixed_mass,
                         ((2. - exc_VE_frac) * fVGdy2(*optres.x) + exc_VE_frac * fVEdy2(*optres.x)) / mass_torsion]],
                       dtype=float)
    w, v = np.linalg.eig(Hessian)
    print(" * First normal mode: {} {} \n   with frequency {} ".format(v[0, 0], v[1, 0], sqrt(w[0]) * 219474.))
    print(" * Second normal mode: {} {} \n   with frequency {} \n".format(v[0, 1], v[1, 1], sqrt(w[1]) * 219474.))
    frequency_torsion.append(sqrt(w[1]) * 219474.)
    frequency_stretching.append(sqrt(w[0]) * 219474.)

if interactiveplots:
    fig, (graph1, graph2) = plt.subplots(2, 1)
    graph1.plot(exc_fraction, frequency_torsion, label="$\omega$ torsion", ls="-", color="r", marker="o")
    graph2.plot(exc_fraction, frequency_stretching, label="$\omega$ stretching", ls="-", color="b", marker="o")
    graph2.set_xlabel("exciton pop at the center")
    graph2.set_ylabel("frequency / cm$^{-1}$")
    graph1.set_ylabel("frequency / cm$^{-1}$")
    graph1.legend()
    graph2.legend()
    plt.show()
    plt.close()
