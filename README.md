# mime

a small Python code to fit a potential energy 
surface to an analytic function defined in a 
sum-of-product form

## License

Copyright (c) 2019, AK Burghardt, Goethe-Universität Frankfurt am Main

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Motivation and concept

The purpose of this code is to derive an analytic 
representation of some data that is dependent on two 
coordinates. This is done by fitting an analytic form
to the data, using least squared. 
In particular, the analytic fitting function is defined
in a sum-over-product form, i.e. two sets of functions
$f_1(x), f_2(x), ..., f_i(x), ...$ and  
$g_1(y), g_2(y), ..., g_i(y), ...$ are defined
for each dimension, and then the fitting function 
is constructed as a linear combination of all the
products: $F(x,y) = \sum_{i,j} a_{i,j} f_i(x) g_i(y)$.

The peculiarity of *mime* is that the lists of functions 
for the two dimensions are not hard-coded, but instead
can be given from input as symbolic expressions.
Thanks to the [sympy](https://www.sympy.org/en/index.html)
library, the symbolic expressions are converted to 
actual python functions. 


## Dependencies

*mime* has been tested with Python3 on an Ubuntu Linux 
distribution. The code depends on a few Python libraries, 
that are listed in the [requirements.txt](requirements.txt)
file. The installation of PyQt5 is suggested for 
visualizing plot to screen in interactive mode with matplotlib,
but is not necessary for the core functions of the code. 


## Installation

For the installation, the testing and the use of **mime**, we 
encourage the use of a Python *virtual environment*

     python3 -m venv mime_test
     source mime_test/bin/activate
     
First, you should get the latest version of the code
cloning the git repository from GitHub

    git clone https://github.com/matbonfanti/mime
    cd mime
    
Then enter in the project root directory and 
install **mime** with the command:

    python setup.py install
    
This command will make the package available and will
install all the necessary requirements.


## Testing

For testing, a unittest class is defined in the 
directory `test`

    cd test
    python test.py
    
This testing unit will evaluate a simple 2D function 
(with additional pseudo-random noise) and then use this
function to fit two different analytic forms. 

The test is passed by checking the values of the coefficients.
The test also produces plots of the original data
(`test_data.pdf`) and of the fitting functions (`plot_fit1.pdf`
and `plot_fit2.pdf`).


## Troubleshooting

In some systems the use of a virtual environment determines 
an **ImportError exception when importing matplotlib**.
This error is determined by the use of the standard Qt5 backend
and can be simply avoided by forcing the choice of 
the matplotlib backend. In practice, we suggest to change
the backend to TkAgg by using the environmental variable MPLBACKEND:

    export MPLBACKEND="TkAgg"

Further details are available in the 
[matplotlib FAQ](https://matplotlib.org/2.1.2/faq/virtualenv_faq.html)


## Examples

The directory `example` contains some realistic use 
of *mime*, to construct an analytic representation of the 
potential energy surfaces for some polymer systems, 
oligothiophene (OT) and para-phenyl vinilene (PPV). 
These functions have been used to study the [electronic 
excitation dynamics in this type of 
system](https://doi.org/10.1063/5.0004510).

