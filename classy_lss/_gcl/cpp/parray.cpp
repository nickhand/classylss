#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "parray.h"

using std::vector;


parray::parray() {
    n[0] = n[1] = n[2] = 0;
}

parray::parray(size_type n0) : vector<double>(n0) {
    n[0] = n0;
    n[1] = n[2] = 1;
}

parray::parray(size_type n0, size_type n1) : vector<double>(n0*n1) {
    n[0] = n0;
    n[1] = n1;
    n[2] = 1;
}

parray::parray(size_type n0, size_type n1, size_type n2) : vector<double>(n0*n1*n2) {
    n[0] = n0;
    n[1] = n1;
    n[2] = n2;
}

parray::parray(size_type n0, const double* v) : vector<double>(v, v+n0) {
    n[0] = n0;
    n[1] = n[2] = 1;
}

parray::parray(size_type n0, size_type n1, const double* v) : vector<double>(v, v+n0*n1) {
    n[0] = n0;
    n[1] = n1;
    n[2] = 1;
}

parray::parray(size_type n0, size_type n1, size_type n2, const double* v) : vector<double>(v, v+n0*n1*n2) {
    n[0] = n0;
    n[1] = n1;
    n[2] = n2;
}

parray::parray(const vector<double>& v) : vector<double>(v) {
    n[0] = (int)size();
    n[1] = n[2] = 1;
}

parray::parray(const parray& v) : vector<double>(v) {
    n[0] = v.n[0];
    n[1] = v.n[1];
    n[2] = v.n[2];
}

parray& parray::operator=(const parray& v) {
    vector<double>::operator=(v);
    n[0] = v.n[0];
    n[1] = v.n[1];
    n[2] = v.n[2];
    return *this;
}

parray::~parray() {
}

parray::parray(size_type n0, size_type n1, size_type n2, double r) : vector<double>(n0*n1*n2, r) {
    n[0] = n0;
    n[1] = n1;
    n[2] = n2;
}

parray parray::zeros(size_type n0) {
    return parray(n0, 1, 1, 0.);
}

parray parray::zeros(size_type n0, size_type n1) {
    return parray(n0, n1, 1, 0.);
}

parray parray::zeros(size_type n0, size_type n1, size_type n2) {
    return parray(n0, n1, n2, 0.);
}

parray parray::ones(size_type n0) {
    return parray(n0, 1, 1, 1.);
}

parray parray::ones(size_type n0, size_type n1) {
    return parray(n0, n1, 1, 1.);
}

parray parray::ones(size_type n0, size_type n1, size_type n2) {
    return parray(n0, n1, n2, 1.);
}

parray::parray(size_type n0, double min, double max, Spacing mode) : vector<double>(n0) {
    n[0] = n0;
    n[1] = n[2] = 1;

    assert(n0 >= 2);
    if(mode == LinearSpacing) {
        for(size_type i = 0; i < n0; i++)
            (*this)[i] = min + i*(max - min)/(n0-1);
    }
    else if(mode == LogarithmicSpacing) {
        assert(min > 0 && max > 0);
        double logmin = log(min);
        double logmax = log(max);
        for(size_type i = 0; i < n0; i++)
            (*this)[i] = exp(logmin + i*(logmax - logmin)/(n0-1));
    }

    /* Guarantee min and max elements are exact */
    (*this)[0] = min;
    (*this)[n0-1] = max;
}

parray parray::linspace(double min, double max, size_type n0) {
    return parray(n0, min, max, LinearSpacing);
}

parray parray::logspace(double min, double max, size_type n0) {
    return parray(n0, min, max, LogarithmicSpacing);
}

void parray::resize(size_type n0) {
    vector<double>::resize(n0);
    n[0] = n0;
    n[1] = n[2] = 1;
}

void parray::resize(size_type n0, size_type n1) {
    vector<double>::resize(n0*n1);
    n[0] = n0;
    n[1] = n1;
    n[2] = 1;
}

void parray::resize(size_type n0, size_type n1, size_type n2) {
    vector<double>::resize(n0*n1*n2);
    n[0] = n0;
    n[1] = n1;
    n[2] = n2;
}

//void parray::transpose(int axis0, int axis1) {
//}

parray parray::slice(int start, int stop, int step) {
    int N = size();
    if (N == 0) return parray();

    if (start < 0) start += N;
    if (stop < 0) stop += N;
    
    // determine the new length
    if (stop > N) stop = N;
    size_type newlen = 1 + (stop-start-1)/step;

    parray r(newlen);
    for (size_type i = 0; i < newlen; i++)
        r[i] = (*this)[start+step*i];
    return r;   
}

parray parray::slice(int start, int stop, int step) const {
    int N = size();
    if (N == 0) return parray();

    if (start < 0) start += N;
    if (stop < 0) stop += N;

    // determine the new length
    if (stop > N) stop = N;
    size_type newlen = 1 + (stop-start-1)/step;

    parray r(newlen);
    for (size_type i = 0; i < newlen; i++)
        r[i] = (*this)[start+step*i];
    return r;   
}

void parray::getshape(size_type* n0, size_type* n1, size_type* n2) const {
    if(n0) *n0 = n[0];
    if(n1) *n1 = n[1];
    if(n2) *n2 = n[2];
}

double parray::min() const {
    size_type N = size();
    if(N == 0)
        return 0.;
    double r = (*this)[0];
    for(size_type i = 1; i < N; i++)
        r = fmin(r, (*this)[i]);
    return r;
}

double parray::max() const {
    size_type N = size();
    if(N == 0)
        return 0.;
    double r = (*this)[0];
    for(size_type i = 1; i < N; i++)
        r = fmax(r, (*this)[i]);
    return r;
}

parray parray::abs() const {
    size_type N = size();
    if(N == 0)
        return parray();
    parray r(N);
    for(size_type i = 0; i < N; i++)
        r[i] = fabs((*this)[i]);
    return r;
}

parray& parray::operator+=(const parray& v) {
    size_type N = size();
    assert(v.size() == N);
    for(size_type i = 0; i < N; i++)
        (*this)[i] += v[i];
    return *this;
}

parray& parray::operator-=(const parray& v) {
    size_type N = size();
    assert(v.size() == N);
    for(size_type i = 0; i < N; i++)
        (*this)[i] -= v[i];
    return *this;
}

parray& parray::operator*=(const parray& v) {
    size_type N = size();
    assert(v.size() == N);
    for(size_type i = 0; i < N; i++)
        (*this)[i] *= v[i];
    return *this;
}

parray& parray::operator/=(const parray& v) {
    size_type N = size();
    assert(v.size() == N);
    for(size_type i = 0; i < N; i++)
        (*this)[i] /= v[i];
    return *this;
}

parray& parray::operator+=(double s) {
    size_type N = size();
    for(size_type i = 0; i < N; i++)
        (*this)[i] += s;
    return *this;
}

parray& parray::operator-=(double s) {
    size_type N = size();
    for(size_type i = 0; i < N; i++)
        (*this)[i] -= s;
    return *this;
}

parray& parray::operator*=(double s) {
    size_type N = size();
    for(size_type i = 0; i < N; i++)
        (*this)[i] *= s;
    return *this;
}

parray& parray::operator/=(double s) {
    size_type N = size();
    for(size_type i = 0; i < N; i++)
        (*this)[i] /= s;
    return *this;
}

parray& parray::pow(double s) {
    size_type N = size();
    for(size_type i = 0; i < N; i++)
        (*this)[i] = std::pow((*this)[i], s);
    return *this;
}

parray operator-(const parray& v) {
    parray::size_type N = v.size();
    parray w(v);
    for(parray::size_type i = 0; i < N; i++)
        w[i] = -w[i];
    return w;
}

parray operator+(const parray& u, const parray& v) {
    assert(u.size() == v.size());
    parray w = u;
    w += v;
    return w;
}

parray operator-(const parray& u, const parray& v) {
    assert(u.size() == v.size());
    parray w = u;
    w -= v;
    return w;
}

parray operator*(const parray& u, const parray& v) {
    assert(u.size() == v.size());
    parray w = u;
    w *= v;
    return w;
}

parray operator/(const parray& u, const parray& v) {
    assert(u.size() == v.size());
    parray w = u;
    w /= v;
    return w;
}

parray operator+(const parray& u, double s) {
    parray w = u;
    w += s;
    return w;
}

parray operator+(double s, const parray& u) {
    parray w = u;
    w += s;
    return w;
}

parray operator-(const parray& u, double s) {
    parray w = u;
    w -= s;
    return w;
}

parray operator-(double s, const parray& u) {
    parray w = u;
    w -= s;
    return -w;
}

parray operator*(const parray& u, double s) {
    parray w = u;
    w *= s;
    return w;
}

parray operator*(double s, const parray& u) {
    parray w = u;
    w *= s;
    return w;
}

parray operator/(const parray& u, double s) {
    parray w = u;
    w /= s;
    return w;
}


parray operator/(double s, const parray& u) {
    parray::size_type N = u.size();
    parray w(N);
    for(parray::size_type i = 0; i < N; i++)
        w[i] = s/u[i];
    return w;
}


double min(const parray& v) {
    return v.min();
}

double max(const parray& v) {
    return v.max();
}

parray abs(const parray& v) {
    return v.abs();
}
