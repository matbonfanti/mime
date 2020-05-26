#!/usr/bin/env python

from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(name='mime',
      version='0.1',
      description='2D function fitting with sum-of-products of symbolic functions',
      long_description=long_description,
      author='Matteo Bonfanti',
      author_email='mat.bonfanti@gmail.com',
      url='https://github.com/matbonfanti/mime',
      packages=['mime'],
      install_requires=["numpy",
                        "matplotlib",
                        "scipy",
                        "setuptools",
                        "sympy",
                        "PyQt5"]
      )
