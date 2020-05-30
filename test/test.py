#! /usr/bin/env python

from math import *
import numpy as np
import mime
import unittest
import create_test_dataset
import os


# =====================================================================================================================

class TestSumOfProductFitting(unittest.TestCase):

    def setUp(self):
        # create test datafile
        create_test_dataset.writedatafile()
        # read and store data from input file
        self.Coord = np.loadtxt(create_test_dataset.DATAFILE, usecols=(0, 1))
        self.Pot = np.loadtxt(create_test_dataset.DATAFILE, usecols=2)

    # ================================================================================================================

    def test_fit1(self):
        # functions for the first degree of freedom
        xfunctions = ["1", "x**1", "x**2", "x**3"]
        # functions for the second degree of freedom
        yfunctions = ["1", "sin(pi*y/180.)", "cos(pi*y/180.)",
                      "sin(2*pi*y/180.)", "cos(2*pi*y/180.)"]

        # define fitting function with all the products between xfunctions and yfunctions
        fittingf = mime.SumOfProduct(xfunctions, yfunctions)
        # do fitting of the data read and stored in self.Coord and self.Pot
        Chisquared = fittingf.fit_data(self.Coord, self.Pot)

        # plot fitting function and residue
        fittingf.plot_function(xLim=create_test_dataset.XLIMITS, nx=30,
                               yLim=create_test_dataset.YLIMITS, ny=30,
                               pdfname="plot_fit1.pdf")

        # The following list has the expected coefficients for the fitting
        targetCoeffs = [2.40933686e-01, 6.06079225e-02, 1.73929492e-01, 1.16803399e-01,
                        -3.64186718e-01, -1.40658799e-01, -5.73230358e-01, -4.80239968e-02,
                        -1.53971857e-01, 9.49311985e-03, 3.86270005e-01, -1.84516952e-01,
                        -8.56936973e-02, 7.73514924e-02, -4.61896573e-04, 5.14052241e-02,
                        -5.14293778e-02, -5.09409547e-02, 2.09858182e-02, -1.17687445e-03]
        # Compare fittingf.Coeffs with targetCoeffs, check if they are equal within numerical precision
        self.assertTrue(np.allclose(targetCoeffs, fittingf.Coeffs, atol=1.e-04))

    # ================================================================================================================

    def test_fit2(self):
        # functions for the first degree of freedom
        xfunctions = ["1", "exp(x)", "exp(2*x)", "exp(3*x)", "exp(-x)"]
        # functions for the second degree of freedom
        yfunctions = ["1", "sinh(pi*y/180.)", "cosh(pi*y/180.)",
                      "sinh(2*pi*y/180.)", "cos(2*pi*y/180.)"]

        # define fitting function with all the products between xfunctions and yfunctions
        fittingf = mime.SumOfProduct(xfunctions, yfunctions)
        # do fitting of the data read and stored in self.Coord and self.Pot
        Chisquared = fittingf.fit_data(self.Coord, self.Pot)

        # plot fitting function and residue
        fittingf.plot_function(xLim=create_test_dataset.XLIMITS, nx=30,
                               yLim=create_test_dataset.YLIMITS, ny=30,
                               pdfname="plot_fit2.pdf")

        # The following list has the expected coefficients for the fitting
        targetCoeffs = [3.68063800e+00, 9.31099956e-01, -3.10498034e+00, 1.96646954e-01,
                        -7.62377398e-01, -2.26436982e+00, -1.19708949e+00, 2.06128180e+00,
                        -6.45191844e-02, 1.69590132e-01, 3.89068690e-01, 2.34145252e-01,
                        -3.28207007e-01, 5.18624042e-03, -4.32760471e-03, -2.51148269e-02,
                        -1.57025699e-02, 2.09746215e-02, -2.88181393e-04, -3.13685857e-04,
                        7.05201105e-01, 1.25803004e-01, -3.99430781e-01, 1.99314539e-02,
                        -5.29416044e-02]
        # Compare fittingf.Coeffs with targetCoeffs, check if they are equal within numerical precision
        self.assertTrue(np.allclose(targetCoeffs, fittingf.Coeffs, atol=1.e-04))

    # ================================================================================================================

    def tearDown(self):
        # clean up memory
        del self.Coord, self.Pot
        # and remove data file
        os.remove("test_data.dat")


# =====================================================================================================================

if __name__ == '__main__':
    unittest.main(verbosity=2)
