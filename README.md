# Bessel_Quadrature_Integrator

Numerical integrals can be hard. Worse than that, numerical integrals can also be very slow. A way to compute some integrals quickly is using Gaussian quadrature rules. In a paper by [Ogata](http://www.kurims.kyoto-u.ac.jp/~okamoto/paper/Publ_RIMS_DE/41-4-40.pdf) in 2005, the author wrote down a quadrature rule for integrals of the form

PICTURE OF MATH

These types of integrals appear very often in cosmology (my field). This package implements that quadrature rule and provides near-optimal implementations of the integral for certain orders of Bessel functions.