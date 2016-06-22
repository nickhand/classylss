#ifndef COMMON_H
#define COMMON_H

#ifndef NULL
#define NULL 0
#endif

// The directory of GCL. This is set to the absolute path to the GCL directory so this is just a failsafe.
#ifndef GCLDIR
#define GCLDIR "." 
#endif

// The directory of GCL. This is set to the absolute path to the data directory so this is just a failsafe.
#ifndef DATADIR
#define DATADIR "." 
#endif

#include <cmath>
#include <cstdio>
#include <string>

/* Forward declarations of generic classes */
class ClassParams;
class ClassEngine;
class parray;
class Cosmology;
class PowerSpectrum;

namespace Common {
  
  
    double SphericalBesselJ0(double x);
    double SphericalBesselJ1(double x);
    double SphericalBesselJ2(double x);
    double SphericalBesselJ3(double x);
    double SphericalBesselJ4(double x);
    double SphericalBesselJ6(double x);
    double SphericalBesselJ8(double x);
    
    /* Convenient print functions that flush the output buffer afterwards */
    //   Print to file
    void write(FILE * stream, const char* format, ...);
    //   Print to stdout
    void info(const char* format, ...);
    //   Print to stdout if VERBOSE is defined.
    void verbose(const char* format, ...);
    //   Print to stdout if DEBUG is defined.
    void debug(const char* format, ...);
    //   Print to stderr.
    void warning(const char* format, ...);
    //   Print to stderr and abort.
    void error(const char* format, ...);
    
    void throw_error(const char *msg, std::string file, int lineno);
    void throw_error(std::string msg, std::string file, int lineno);


    /***** Math routines *****/

    /* Small integer powers */
    static inline double pow2(double x) { return x*x; }
    static inline double pow3(double x) { return x*x*x; }
    static inline double pow4(double x) { return pow2(pow2(x)); }
    static inline double pow5(double x) { return x*pow4(x); }
    static inline double pow6(double x) { return pow3(pow2(x)); }
    static inline double pow7(double x) { return x*pow6(x); }
    static inline double pow8(double x) { return pow2(pow4(x)); }
    static inline double pow9(double x) { return x*pow8(x); }
    static inline double pow10(double x) { return pow2(pow5(x)); }
}

/* Physical constants in cgs units */
namespace Constants {
    
    // set up cgs units
    const double cm     = 1.;
    const double gram   = 1.;
    const double second = 1.;
    const double erg    = 1.;
    const double kelvin = 1.;
    const double radian = 1.;

    // prefixes
    const double giga  = 1e9;
    const double mega  = 1e6;
    const double kilo  = 1e3;
    const double centi = 1e-2;
    const double milli = 1e-3;
    const double micro = 1e-6;
    const double nano  = 1e-9;
    const double pico  = 1e-12;

    // fundamental constants
    const double h_planck = 6.62606957e-27;     /* planck constant */
    const double h_bar    = h_planck/(2.*M_PI); /* reduced planck constant */
    const double c_light  = 2.99792458e10;      /* speed of light */
    const double k_b      = 1.3806488e-16;      /* Boltzmann constant */
    const double m_p      = 1.672621777e-24;    /* proton mass */
    const double m_e      = 9.10938291e-28;     /* electron mass */
    const double q_e      = 4.80320425e-10;     /* e.s.u. */
    const double G        = 6.67384e-8;         /* Newton's constant */
    const double eV       = 1.60217657e-12;     /* 1 eV in ergs */
    const double N_a      = 6.02214129e23;      /* Avogadro's constant */
    const double sigma_sb = 5.670373e-5;        /* Stefan-Boltzmann constant */
    const double a_rad    = 4*sigma_sb/c_light; /* radiation constant */
    const double sigma_T  = 6.652458734e-25;    /* Thomson cross section */
    const double T_cmb    = 2.72528;            /* temperature of the CMB */
    const double H_0      = 100.;               /* in units of h km/s/Mpc */
    const double a_0      = 5.2917721092e-9;    /* Bohr radius */

    // length conversion factors
    const double km       = 1e5 * cm;  
    const double meter    = 1e2 * cm;  
    const double inch     = 2.54 * cm;  
    const double mm       = 1e-1 * cm;  
    const double micron   = 1e-4 * cm;  
    const double angstrom = 1e-8 * cm;
    const double jansky   = 1e-23;       /* in erg/s/cm/cm/Hz */
    const double barn     = 1e-24*cm*cm;

    // energy and power conversion factors
    const double joule   = 1e7 * erg;
    const double watt    = 1e7 * erg / second;
    const double rydberg = 10973731.568539 * 1./meter * h_planck*c_light; /* in ergs */

    // astronomical constants
    const double au        = 149597870700.*meter;   /* the astronomical unit */
    const double degree    = M_PI/180.;
    const double arcminute = degree/60.;
    const double arcsecond = arcminute/60.;         /* radians */
    const double parsec    = au/arcsecond;          /* cm */
    const double minute    = 60.;                    /* seconds */
    const double hour      = 60.*60.;              /* seconds */
    const double day       = 8.64e4;                 /* seconds */
    const double year      = 365.2425 * day;         /* seconds */
    const double Mpc       = mega*parsec;
    const double lyr       = c_light*year;           /* light year */

    const double L_sun = 3.826e33;                 /* erg/s */
    const double M_sun = 1.9891e33;                /* g */
    const double R_sun = 6.9598e10;                /* cm */
    const double T_sun = 5770;                     /* kelvin */

    const double M_earth   = 5.976e27;              /* Earth mass in g */
    const double R_earth   = 6371 * km;             /* Earth's equatorial radius in cm */
    const double M_jupiter = 1898.8e27;             /* Jupiter mass in g */
    const double R_jupiter = 70850 * km;            /* Jupiter's equatorial radius in cm */

}

#endif // COMMON_H
