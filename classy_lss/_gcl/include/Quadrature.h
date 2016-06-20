/**
 * @file Quadrature.h
 * @author Jordan Carlson (jwgcarlson@gmail.com)
 * @brief Templated routines for numerical quadrature. */

#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <cmath>
#include "parray.h"

double TrapzIntegrate(const parray& x, const parray& y);
double SimpsIntegrate(const parray& x, const parray& y);
double basic_simps(const parray& y, int start, int stop, const parray& x);

/* General notes on integration routines:
 *  - The Function f can be any callable object, i.e. any object that defines
 *    the method  double operator()(double x).  This becomes especially
 *    powerful when combined with the Boost/TR1 bind() method.  Of course you
 *    can also use regular C-style functions of signature  double f(double x).
 *  - epsrel is the desired relative error, epsabs is the desired absolute
 *    error.  The adaptive integration routines continue until _either_
 *      relative error < epsrel or absolute error < epsabs.
 *  - abserr is the computed absolute error.
 *  - neval is the number of times the integrand was evaluated.
 *  - Substitution mixins make it easy to change the integration variable
 *    (usually to improve convergence) without having to redefine the function
 *    itself. */


/**
 * \brief Integrate \a f(x) from \a a to \a b.
 *
 * Compute the definite integral $\int_a^b f(x) dx$. 
 *
 * Uses an adaptive algorithm with a 15-point Gauss-Kronrod rule.  Based on
 * the implementation of gsl_integration_qag in the GNU Scientific Library. */
template<typename Function>
double Integrate(Function f, double a, double b, double epsrel = 1e-5, double epsabs = 1e-10, double* abserr = 0, int* neval = 0);

/**
 * \brief Integrate \a f(x) from \a a to \a b using the given substition rule.
 *
 * Compute the definite integral $\int_a^b f(x) dx$ by making the change of
 * variables x -> u, i.e. compute $\int_{u(a)}^{u(b)} f(x(u)) (dx/du) du$.
 * This can be used to speed up the convergence of an integral by choosing an
 * substitution for which the integrand appears smooth. */
template<typename Sub, typename Function>
double Integrate(Function f, double a, double b, double epsrel = 1e-5, double epsabs = 1e-10, double* abserr = 0, int* neval = 0, Sub sub = Sub());

struct NoSub { \
    double x(double u) { return u; } \
    double u(double x) { return x; } \
    double dxdu(double) { return 1; } \
};

struct ExpSub { \
    double x(double u) { return exp(u); } \
    double u(double x) { return log(x); } \
    double dxdu(double u) { return exp(u); } \
};

struct InverseSub { \
    double x(double u) { return 1/u; } \
    double u(double x) { return 1/x; } \
    double dxdu(double u) { return -1/(u*u); } \
};


/**
 * \brief Compute an n-dimensional definite integral.
 *
 * Compute $\int f(\vec{x}) d^nx$ using the adaptive subregion algorithm of
 * Genz and Malik.  The integration region is the rectangle defined by the
 * parrays a and b, i.e.
 *   V = [a_1,b_1] x [a_2,b_2] x ... x [a_n,b_n].
 * For n <= 6, the function f should have the signature
 *   double f(double x1, ..., double xn);
 * For n > 6 the signature should be
 *   double f(double* x);
 * The Function template parameter can be determined implicitly from f, so you
 * need not specify it.  However, the dimensionality n _is_ required, so you
 * must write e.g.
 *   double s = Integrate<n>(f, a, b);
 *
 * This implementation is based on the Fortran routine ADAPT from Genz & Malik,
 * J. Comp. & Appl. Math., 6, 295-302, 1980. */
template<int n, typename Function>
double Integrate(Function f, double* a, double* b, double epsrel = 1e-5, double epsabs = 1e-10, double* abserr = 0, int* neval = 0);


#include "Quadrature.inl"

#endif // QUADRATURE_H
