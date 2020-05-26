#! /usr/bin/env python

import math
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from write_mctdh_opfile import *
import mime

interactiveplots = False
pdfplots = True

# #################################################
#              2D V_G potential Fit               #
# #################################################

print("\n *** 2D V_G potential fit ***\n")

# functions for the first degree of freedom (write functions of x)
xFunctLabels = ["1", "x**1", "x**2", "x**3"]
xMCTDHLabels = ["1", "q^1", "q^2", "q^3"]

# functions for the second degree of freedom (write functions of y)
yFunctLabels = ["1", "cos(y)", "cos(2*y)", "cos(3*y)", "cos(4*y)"]
yMCTDHLabels = ["1", "cos[1.0,0.0]", "cos[2.0,0.0]", "cos[3.0,0.0]", "cos[4.0,0.0]"]

# define fitting function
VG_2D = mime.sumofproduct(xFunctLabels, yFunctLabels, xBasisFormat=xMCTDHLabels, yBasisFormat=yMCTDHLabels)

# file to read the input data
inputfile = "ot-fitdata.dat"
# read and store data from input file
Coord = np.loadtxt(inputfile, usecols=(3, 0))
Pot = np.loadtxt(inputfile, usecols=4)
# convert angles from deg to rad
Coord = np.array([[x, (theta - 180.0) * math.pi / 180.0] for x, theta in Coord])

# do fitting and get optimal coefficients
Chisquared = VG_2D.fit_data(Coord, Pot, FixZero=False)
NDoF = len(Pot) - len(xFunctLabels) * len(yFunctLabels)
print(" Fitting done... regression standard error is {} meV".format(math.sqrt(Chisquared / NDoF) * 27.21142))

print("\n{:>15s} | {:60s}".format("coeff.s", "basis functions"))
for c, f in zip(VG_2D.get_coeffs(), VG_2D.get_basis_funcs()):
    print("{:15.6e} | {:60s}".format(c, str(f)))

# visualize data and model 
if interactiveplots:
    VG_2D.show_1D_modelplusdata(Coord, Pot, yspacing=180.0 / math.pi)

# plot graphs
if pdfplots:
    Name = "VG_2D"
    print("\n Writing potential to file {}".format('"' + Name + '_func.pdf"'))
    VG_2D.plot_function(xLim=(-0.3, +0.3), yLim=(-0.5 * math.pi, 0.5 * math.pi), nx=27, ny=37,
                        pdfname=Name + "_func.pdf", scale_y=180.0 / math.pi, scale_v=27.21142)
    print(" Writing fitting residues to file {}\n".format('"' + Name + '_residue.pdf"'))
    VG_2D.plot_residue(Coord, Pot, nx=27, ny=19, pdfname=Name + "_residue.pdf", scale_y=180.0 / math.pi,
                       scale_v=27.21142)

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
inputfile = "ot-fitdata.dat"
# read and store data from input file
Coord = np.loadtxt(inputfile, usecols=(3, 0))
Pot = np.loadtxt(inputfile, usecols=5)
# convert angles from deg to rad
Coord = np.array([[x, (theta - 180.0) * math.pi / 180.0] for x, theta in Coord])

# do fitting and get optimal coefficients
Chisquared = VE_2D.fit_data(Coord, Pot, FixZero=False)
NDoF = len(Pot) - len(xFunctLabels) * len(yFunctLabels)
print(" Fitting done... regression standard error is {} meV".format(math.sqrt(Chisquared / NDoF) * 27.21142))

print("\n{:>15s} | {:60s}".format("coeff.s", "basis functions"))
for c, f in zip(VE_2D.get_coeffs(), VE_2D.get_basis_funcs()):
    print("{:15.6e} | {:60s}".format(c, str(f)))

# visualize data and model 
if interactiveplots:
    VE_2D.show_1D_modelplusdata(Coord, Pot, yspacing=180.0 / math.pi)

# plot graphs
if pdfplots:
    Name = "VE_2D"
    print("\n Writing potential to file {}".format('"' + Name + '_func.pdf"'))
    VE_2D.plot_function(xLim=(-0.3, +0.3), yLim=(-0.5 * math.pi, 0.5 * math.pi), nx=27, ny=37,
                        pdfname=Name + "_func.pdf", scale_y=180.0 / math.pi, scale_v=27.21142)
    print(" Writing fitting residues to file {}\n".format('"' + Name + '_residue.pdf"'))
    VE_2D.plot_residue(Coord, Pot, nx=27, ny=19, pdfname=Name + "_residue.pdf", scale_y=180.0 / math.pi,
                       scale_v=27.21142)

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
inputfile = "ot-fitdata.dat"
# read and store data from input file
Coord = np.loadtxt(inputfile, usecols=(3, 0))
Pot = np.loadtxt(inputfile, usecols=6)
# convert angles from deg to rad
Coord = np.array([[x, (theta - 180.0) * math.pi / 180.0] for x, theta in Coord])

# do fitting and get optimal coefficients
Chisquared = w_2D.fit_data(Coord, Pot)
NDoF = len(Pot) - len(xFunctLabels) * len(yFunctLabels)
print(" Fitting done... regression standard error is {} meV".format(math.sqrt(Chisquared / NDoF) * 27.21142))

print("\n{:>15s} | {:60s}".format("coeff.s", "basis functions"))
for c, f in zip(w_2D.get_coeffs(), w_2D.get_basis_funcs()):
    print("{:15.6e} | {:60s}".format(c, str(f)))

# visualize data and model 
if interactiveplots:
    w_2D.show_1D_modelplusdata(Coord, Pot, yspacing=180.0 / math.pi)

# plot graphs
if pdfplots:
    Name = "w_2D"
    print("\n Writing potential to file {}".format('"' + Name + '_func.pdf"'))
    w_2D.plot_function(xLim=(-0.3, +0.3), yLim=(-0.5 * math.pi, 0.5 * math.pi), nx=27, ny=37,
                       pdfname=Name + "_func.pdf", scale_y=180.0 / math.pi, scale_v=27.21142)
    print(" Writing fitting residues to file {}\n".format('"' + Name + '_residue.pdf"'))
    w_2D.plot_residue(Coord, Pot, nx=27, ny=19, pdfname=Name + "_residue.pdf", scale_y=180.0 / math.pi,
                      scale_v=27.21142)

# #################################################o

print("\n *** 1D ring breathing ground and exc states potential fit ***\n")


def MorsePotential1(Crd, D, alpha):
    return D * (np.exp(-alpha * Crd) - 1.0) ** 2


def MorsePotential2(Crd, D, alpha, x0, E0):
    return D * (np.exp(-alpha * (Crd - x0)) - 1.0) ** 2 + E0


# file to read the input data
inputfile = "ot-ring_VG_VE.dat"
# read and store data from input file
Coord = np.loadtxt(inputfile, usecols=0)
PotGS = np.loadtxt(inputfile, usecols=1)
PotES = np.loadtxt(inputfile, usecols=2)

RingGS_Morse, _ = curve_fit(MorsePotential1, Coord, PotGS, p0=(1.0, 1.0))
RingES_Morse, _ = curve_fit(MorsePotential2, Coord, PotES, p0=(1.0, 1.0, 0.0, 0.0))

if interactiveplots:
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

print("\n *** MCTDH operator ***\n")

NameOperator = "20-mer.op"
print(" Writing operator file to file {} ".format(NameOperator))

write_mctdh_opfile(20, 20, VG_2D, VE_2D, w_2D, RingGS_Morse, RingES_Morse)
