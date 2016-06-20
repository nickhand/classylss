#include <cassert>
#include <cstdio>
#include <ctime>
#include <functional>

#include "Quadrature.h"
#include "PowerSpectrum.h"

using std::bind;
using std::cref;
using namespace std::placeholders;
using namespace Common;

PowerSpectrum::PowerSpectrum() {}

PowerSpectrum::~PowerSpectrum() {}

parray PowerSpectrum::EvaluateMany(const parray& k) const {
    int n = (int)k.size();
    parray pk(n);
    #pragma omp parallel for
    for(int i = 0; i < n; i++)
        pk[i] = Evaluate(k[i]);
    return pk;
}

/* Top-hat window function */
static double W(double x) {
    if(x < 1e-5)
        return 1 - (1./30.)*x*x;
    else
        return 3/pow3(x) * (sin(x) - x*cos(x));
}

static double f(const PowerSpectrum& P, double R, double k) {
    return k*k/(2*M_PI*M_PI) * P(k) * pow2(W(k*R));
}
double PowerSpectrum::Sigma(double R) const {
    double sigma2 = Integrate<ExpSub>(bind(f, cref(*this), R, _1), 1e-5, 1e2, 1e-5, 1e-12);
    return sqrt(sigma2);
}

parray PowerSpectrum::Sigma(const parray& R) const {
    int n = (int)R.size();
    parray sig(n);
    #pragma omp parallel for
    for(int i = 0; i < n; i++)
        sig[i] = Sigma(R[i]);
    return sig;
}

static double g(const PowerSpectrum& P, double k) {
    return P(k);
}
double PowerSpectrum::VelocityDispersion() const {
    return 1/(6*M_PI*M_PI) * Integrate<ExpSub>(bind(g, cref(*this), _1), 1e-5, 1e2, 1e-5, 1e-12);
}

double PowerSpectrum::VelocityDispersion(double k, double factor ) const {
    return 1/(6*M_PI*M_PI) * Integrate<ExpSub>(bind(g, cref(*this), _1), 1e-5, factor*k, 1e-5, 1e-12);
}

parray PowerSpectrum::VelocityDispersion(const parray& k, double factor) const {
    int n = (int)k.size();
    parray sigmasq(n);
    #pragma omp parallel for
    for(int i = 0; i < n; i++)
        sigmasq[i] = VelocityDispersion(k[i], factor);
    return sigmasq;
}

double PowerSpectrum::NonlinearScale() const {
    return 1/sqrt(VelocityDispersion());
}


