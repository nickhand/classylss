#ifndef PARRAY_H
#define PARRAY_H

#include <vector>

#include "Common.h"

/* A simple multi-dimensional (up to 3) array, based on Numpy's array type.
 * The multi-dimensional array is represented in memory as a flat array
 * inherited from the STL vector class.  For 1-D arrays it may be used as
 * a plain old std::vector<double>. For 2-D and 3-D arrays, the storage order
 * is row major. */
class parray : public std::vector<double> {

public:
    /* Construct an empty (0-dimensional) array. */
    parray();

    /* Construct an array with arbitrary initial values. */
    explicit parray(size_type n0);
    parray(size_type n0, size_type n1);
    parray(size_type n0, size_type n1, size_type n2);

    /* Construct a 1-D array of length n, with all elements set to r. */
    // array(size_type n, double r);

    /* Construct an array with elements initialized from the C-array v.  For
     * 2-D arrays, the data elements are presumed to be laid out in memory as
     *   a_{ij} = v[n1*i + j] ,
     * i.e. in row-major order.  For 3-D arrays the memory layout is
     *   a_{ijk} = v[n2*(n1*i + j) + k] .  */
    parray(size_type n0, const double* v);
    parray(size_type n0, size_type n1, const double* v);
    parray(size_type n0, size_type n1, size_type n2, const double* v);

    /* Construct a 1-D array from an STL vector */
    parray(const std::vector<double>& v);

    /* Copy constructor. */
    parray(const parray& v);

    /* Assignment.  This array is resized, and the elements of v are copied. */
    parray& operator=(const parray& v);

    /* Destructor. */
    ~parray();

    /* Construct an array with values initialized to zero. */
    static parray zeros(size_type n0);
    static parray zeros(size_type n0, size_type n1);
    static parray zeros(size_type n0, size_type n1, size_type n2);

    /* Construct an array with values initialized to one. */
    static parray ones(size_type n0);
    static parray ones(size_type n0, size_type n1);
    static parray ones(size_type n0, size_type n1, size_type n2);

    /* Constuct a 1-D array of length n with linearly or logarithmically spaced
     * values between min and max.  Note that if min is greater than max, the
     * values will be arranged in descending order.  Also note that for
     * logspace(), both min and max must be strictly positive. */
    static parray linspace(double min, double max, size_type n = 50);
    static parray logspace(double min, double max, size_type n = 50);

    /* Resize this array. */
    void resize(size_type n0);
    void resize(size_type n0, size_type n1);
    void resize(size_type n0, size_type n1, size_type n2);

    /* Transpose the current array over two axes. */
    // void transpose(int axis0 = 0, int axis1 = 1);

    void getshape(size_type* n0, size_type* n1, size_type* n2) const;

    double min() const;
    double max() const;
    parray abs() const;

    /* Standard flat accessors. (allows for negative indexing!) */
    const double& operator[](int i) const { if (i<0) i += size(); return std::vector<double>::operator[](i); }
    double& operator[](int i) { if (i<0) i += size(); return std::vector<double>::operator[](i); }
    
    /* Slicing operators */
    parray slice(int start, int stop, int step=1) const;
    parray slice(int start, int stop, int step=1);

    /* Accessors for 1-D arrays */
    const double& operator()(int i0) const { return (*this)[i0]; }
    double& operator()(int i0) { return (*this)[i0]; }

    /* Accessors for 2-D arrays */
    const double& operator()(int i0, int i1) const { return (*this)[i0*n[1] + i1]; }
    double& operator()(int i0, int i1) { return (*this)[i0*n[1] + i1]; }

    /* Accessors for 3-D arrays */
    const double& operator()(int i0, int i1, int i2) const { return (*this)[(i0*n[1] + i1)*n[2] + i2]; }
    double& operator()(int i0, int i1, int i2) { return (*this)[(i0*n[1] + i1)*n[2] + i2]; }

    /* Automatic cast to C-array.  Note that the returned pointer is only valid
     * as long as this array is not resized. */
    operator double*() { return data(); }
    operator const double*() const { return data(); }

    /* Arithmetic operations.  Note that these methods do not check the shape
     * of the array v, only that its total size equals the total size of this
     * array. */
    parray& operator+=(const parray& v);  // element-wise addition
    parray& operator-=(const parray& v);  // element-wise subtraction
    parray& operator*=(const parray& v);  // element-wise multiplication
    parray& operator/=(const parray& v);  // element-wise division
    parray& operator+=(double s);        // add a constant from each element
    parray& operator-=(double s);        // subtract a constant from each element
    parray& operator*=(double s);        // multiplication by a constant
    parray& operator/=(double s);        // division by a constant
    parray& pow(double s);

protected:
    size_type n[3];     // shape of array

    enum Spacing {
        LinearSpacing,
        LogarithmicSpacing
    };
    parray(size_type n0, double min, double max, Spacing mode);

    parray(size_type n0, size_type n1, size_type n2, double r);
};

parray operator-(const parray& v);    // unary negation
parray operator+(const parray& u, const parray& v);
parray operator-(const parray& u, const parray& v);
parray operator*(const parray& u, const parray& v);
parray operator/(const parray& u, const parray& v);
parray operator+(const parray& u, double s);
parray operator+(double s, const parray& u);
parray operator-(const parray& u, double s);
parray operator-(double s, const parray& u);
parray operator*(const parray& u, double s);
parray operator*(double s, const parray& u);
parray operator/(const parray& u, double s);
parray operator/(double s, const parray& u);

double min(const parray& v);
double max(const parray& v);
parray abs(const parray& v);

#endif // PARRAY_H
