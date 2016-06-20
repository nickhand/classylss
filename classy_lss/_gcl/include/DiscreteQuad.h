
#ifndef DISCRETEQUAD_H
#define DISCRETEQUAD_H

#include <cmath>
#include "parray.h"


/* mirrors scipy.integrate.trapz function */
double TrapzIntegrate(const parray& x, const parray& y);

/* mirrors scipy.integrate.simps function, with `even` keyword equal to `avg` */
double SimpsIntegrate(const parray& x, const parray& y);
double basic_simps(const parray& y, int start, int stop, const parray& x);

/**
 * \brief Integrate f(x) based on its values at n uniformly spaced points.
 *
 * Compute a weighted summation approximating the definite integral $\int f(x)
 * dx$ for a function f(x) that has already been evaluated at n discrete points
 * with uniform spacing h.  Uses Simpson's rule for n odd, Hollingsworth and
 * Hunter's 3rd-order formula for n even. */
double DiscreteIntegrate(int n, double* f, double h = 1);

#endif
