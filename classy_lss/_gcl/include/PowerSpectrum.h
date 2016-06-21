#ifndef POWER_SPECTRUM_H
#define POWER_SPECTRUM_H

/* Conventions for power spectra:
 *   - k is measured in h/Mpc
 *   - P(k) is measured in (Mpc/h)^3
 *   - My Fourier conventions are
 *       \tilde f(k) = \int dx e^{-ikx} f(x)
 *       f(x) = \int \frac{dk}{(2\pi)^3} e^{ikx} \tilde f(k)
 *     which means P(k) is defined by
 *      (2pi)^3 \delta_D(k+k') P(k) = < \delta(k) \delta(k') >
*/

#include "Common.h"
#include "parray.h"


/**********************************************************************
 * PowerSpectrum
 *
 * Base class for double-space power spectra.
 *********************************************************************/

class PowerSpectrum {
public:
    PowerSpectrum();
    virtual ~PowerSpectrum();

    virtual const Cosmology& GetCosmology() const = 0;

    /* Evaluate power spectrum at a single k */
    virtual double Evaluate(double k) const = 0;
    double operator()(double k) const { return Evaluate(k); }

    /* Evaluate power spectrum at many k values (parallelized for speed) */
    virtual parray EvaluateMany(const parray& k) const;
    parray operator()(const parray& k) const { return EvaluateMany(k); }

    /* Calculate the variance with top-hat smoothing radius R Mpc/h */
    virtual double Sigma(double R) const;
    virtual parray Sigma(const parray& R) const;

    /* Calculate the 1-D velocity dispersion $\sigma_v^2 = \frac{1}{3} \int_0^\infty P(q) ~dq$ */
    virtual double VelocityDispersion() const;
    
    /* Calculate the 1-D velocity dispersion $\sigma_v^2(k) = \frac{1}{3} \int_0^{factor*k} P(q) ~dq$ */
    virtual double VelocityDispersion(double k, double factor = 0.5) const;
    
    /* compute \sigma_v^2(k) for several k */
    virtual parray VelocityDispersion(const parray& k, double factor = 0.5) const;

    /* Calculate the non-linear wavenumber scale $k_\text{nl} = 1/\sigma_v$ */
    virtual double NonlinearScale() const;
    
    virtual const double& GetRedshift() const = 0;
    
    /* Zel'dovich approx integral X(k) = 1/2pi^2 \int dq P(q) (-2 * j_1(kq) / kq) */
    virtual double X_Zel(double k) const;
    
    /* Compute X(k) for several k*/
    virtual parray X_Zel(const parray& k) const;
    
    /* Zel'dovich approx integral Y(k) = 1/2pi^2  \int dq P(q) (-2 * j_0(kq) + 6 j_1(kq)/kq) */
    virtual double Y_Zel(double k) const;
    
    /* Compute Y(k) for several k */
    virtual parray Y_Zel(const parray& k) const;
    
    /* Zel'dovich integral Q3(k) = 1/(10pi^2) k^4 \int dq P(q)^2/q^2 */
    virtual double Q3_Zel(double k) const;
    
    /* Compute Q3(k) for several k */
    virtual parray Q3_Zel(const parray& k) const;
};


#endif /* POWER_SPECTRUM_H */
