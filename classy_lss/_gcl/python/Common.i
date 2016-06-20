%{
#include "Common.h"
%}

/* Type definitions */
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

/* Physical constants in SI units */
%nodefaultctor Constants;
%nodefaultdtor Constants;
struct Constants {
    // set up cgs units
    static const double cm     = 1.;
    static const double gram   = 1.;
    static const double second = 1.;
    static const double erg    = 1.;
    static const double kelvin = 1.;
    static const double radian = 1.;

    // prefixes
    static const double giga  = 1e9;
    static const double mega  = 1e6;
    static const double kilo  = 1e3;
    static const double centi = 1e-2;
    static const double milli = 1e-3;
    static const double micro = 1e-6;
    static const double nano  = 1e-9;
    static const double pico  = 1e-12;

    // fundamental constants
    static const double h_planck = 6.62606957e-27;     /* planck constant */
    static const double h_bar    = h_planck/(2.*M_PI); /* reduced planck constant */
    static const double c_light  = 2.99792458e10;      /* speed of light */
    static const double k_b      = 1.3806488e-16;      /* Boltzmann constant */
    static const double m_p      = 1.672621777e-24;    /* proton mass */
    static const double m_e      = 9.10938291e-28;     /* electron mass */
    static const double q_e      = 4.80320425e-10;     /* e.s.u. */
    static const double G        = 6.67384e-8;         /* Newton's constant */
    static const double eV       = 1.60217657e-12;     /* 1 eV in ergs */
    static const double N_a      = 6.02214129e23;      /* Avogadro's constant */
    static const double sigma_sb = 5.670373e-5;        /* Stefan-Boltzmann constant */
    static const double a_rad    = 4*sigma_sb/c_light; /* radiation constant */
    static const double sigma_T  = 6.652458734e-25;    /* Thomson cross section */
    static const double T_cmb    = 2.72528;            /* temperature of the CMB */
    static const double H_0      = 100.;               /* in units of h km/s/Mpc */
    static const double a_0      = 5.2917721092e-9;    /* Bohr radius */

    // length conversion factors
    static const double km       = 1e5 * cm;  
    static const double meter    = 1e2 * cm;  
    static const double inch     = 2.54 * cm;  
    static const double mm       = 1e-1 * cm;  
    static const double micron   = 1e-4 * cm;  
    static const double angstrom = 1e-8 * cm;
    static const double jansky   = 1e-23;       /* in erg/s/cm/cm/Hz */
    static const double barn     = 1e-24*cm*cm;

    // energy and power conversion factors
    static const double joule   = 1e7 * erg;
    static const double watt    = 1e7 * erg / second;
    static const double rydberg = 10973731.568539 * 1./meter * h_planck*c_light; /* in ergs */

    // astronomical constants
    static const double au        = 149597870700.*meter;   /* the astronomical unit */
    static const double degree    = M_PI/180.;
    static const double arcminute = degree/60.;
    static const double arcsecond = arcminute/60.;         /* radians */
    static const double parsec    = au/arcsecond;          /* cm */
    static const double minute    = 60.;                    /* seconds */
    static const double hour      = 60.*60.;              /* seconds */
    static const double day       = 8.64e4;                 /* seconds */
    static const double year      = 365.2425 * day;         /* seconds */
    static const double Mpc       = mega*parsec;
    static const double lyr       = c_light*year;           /* light year */

    static const double L_sun = 3.826e33;                 /* erg/s */
    static const double M_sun = 1.9891e33;                /* g */
    static const double R_sun = 6.9598e10;                /* cm */
    static const double T_sun = 5770;                     /* kelvin */

    static const double M_earth   = 5.976e27;              /* Earth mass in g */
    static const double R_earth   = 6371 * km;             /* Earth's equatorial radius in cm */
    static const double M_jupiter = 1898.8e27;             /* Jupiter mass in g */
    static const double R_jupiter = 70850 * km;            /* Jupiter's equatorial radius in cm */
};


/* Allow any float-able type to be treated as a double */
%typemap(in) double {
    PyObject* floatobj = PyNumber_Float($input);
    $1 = PyFloat_AsDouble(floatobj);
    Py_DECREF(floatobj);
}
%typemap(typecheck) double {
    // check for numpy array first, so we don't raise a TypeError later on
    if (($input) && PyArray_Check((PyArrayObject *)$input))
        $1 = 0;
    else {
        PyObject* floatobj = PyNumber_Float($input);
        if(floatobj != NULL) {
            $1 = 1;
            Py_DECREF(floatobj);
        }
        else 
            $1 = 0;
    }
}
