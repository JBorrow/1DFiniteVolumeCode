1D Finite Volume Code
=====================

This code was produced as part of the [SAMCRHSS2019](http://www-star.st-and.ac.uk/samcss/)
at the University of St. Andrews in August 2019. This is a 1D finite element code, produced
using the information from Bert Vandenbroucke's lectures, and includes both the
first and second order variants. It is written in pure C99, and uses the exact Riemann
solver from the [SWIFT](http://swiftsim.com) code.

To build the code and run it, all you need is a C compiler and python with the matplotlib
and numpy libraries. You can do this with
```
make all
```
which will produce two plots, one showing the Sod Shock in 1D with the first order code, and
the second showing it with the second order code.

Hopefully this code is at least somewhat well strucutured and can be used by others to learn
more about these methods.