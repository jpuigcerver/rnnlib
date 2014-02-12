RNNLIB
======

Disclaimer
----------

I found difficult going through the Grave's RNNLIB code. This is an attempt
of making it more readable, following the Google C++ Style guide
(http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml) and adding
some comments.

The original code was downloaded from the RNNLIB official webpage
(date 02/14/14):

- http://sourceforge.net/p/rnnl/wiki/Home/

For cutting edge versions, refer to that webpage. I will try to keep this
updated, but I can't make any promise.

The original project was released under a GPLv3 license, and the same applies
here. All credit and bug reports should go to Alex Graves.


ORIGINAL README
===============

RNNLIB is written in C++ and should compile on any platform. However it is
currently only tested for Linux and OSX.

Building it requires the following:

* A modern C++ compiler (e.g. gcc 3.0 or higher)
* GNU Libtool
* GNU automake version 1.9 or higher
* the NetCDF scientific data library
* Boost C++ Libraries version 1.36 or higher (headers only, no compilation
needed)

In addition, the following python packages are needed for the auxiliary scripts
in the 'utils' directory:

* SciPy
* PyLab
* PIL

And this package is needed to create and manipulate netcdf data files with
python, and to run the experiments in the 'examples' directory:

* ScientificPython (NOT Scipy)

To build RNNLIB, first download the source from here, then enter the root
directory and type

$ ./configure
$ make

This should create the binary file 'rnnlib' in the 'bin' directory.
Note that on most linux systems the default installation directory for the
Boost headers is '/usr/local/include/boost-VERSION_NUMBER' which is not on the
standard include path.
In this case type

$ CXXFLAGS=-I/usr/local/include/boost-VERSION_NUMBER/ ./configure
$ make

If you wish to install the binary type:

$ make install

By default this will use '/usr' as the installation root (for which you will
usually need administrator privileges).
You can change the install path with the —prefix option of the configure script
(use ./configure —help for other options)

It is recommended that you add the directory containing the 'rnnlib' binary to
your path, as otherwise the tools in the 'utilities' directory will not work.

Project files are provided for the following integrated development
environments in the 'ide' directory:

* kdevelop (KDE, linux)
* xcode (OSX)
