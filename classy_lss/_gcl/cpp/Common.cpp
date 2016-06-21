#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <stdexcept>

#include "Common.h"

double Common::SphericalBesselJ0(double x) {
    if(fabs(x) < 1e-4)
        return 1 - Common::pow2(x)/6 + Common::pow4(x)/120 - Common::pow6(x)/5040;
    else
        return sin(x)/x;
}

double Common::SphericalBesselJ1(double x) {
    if(fabs(x) < 1e-4)
        return x/3 - Common::pow3(x)/30 + Common::pow5(x)/840 - Common::pow7(x)/45360;
    else
        return (sin(x) - x*cos(x))/Common::pow2(x);
}

double Common::SphericalBesselJ2(double x) {
    if(fabs(x) < 1e-4)
        return Common::pow2(x)/15 - Common::pow4(x)/210 + Common::pow6(x)/7560;
    else
        return ((3 - Common::pow2(x))*sin(x) - 3*x*cos(x))/Common::pow3(x);
}

double Common::SphericalBesselJ3(double x) {
    if(fabs(x) < 1e-4)
        return Common::pow3(x)/105 - Common::pow5(x)/1890 + Common::pow7(x)/83160;
    else
        return ((15 - 6*Common::pow2(x))*sin(x) - (15*x - Common::pow3(x))*cos(x))/Common::pow4(x);
}

double Common::SphericalBesselJ4(double x) {
    if(fabs(x) < 1e-4)
        return Common::pow4(x)/945 - Common::pow6(x)/20790;
    else
        return 5.*(2.*Common::pow2(x) - 21.)*cos(x)/Common::pow4(x) + (Common::pow4(x) - 45.*Common::pow2(x) + 105.)*sin(x)/Common::pow5(x);
}

double Common::SphericalBesselJ6(double x) {
    if(fabs(x) < 1e-4)
        return Common::pow6(x)/135135 - Common::pow8(x)/4054050;
    else
        return ((-Common::pow6(x) + 210*Common::pow4(x) - 4725*Common::pow2(x) + 10395)*sin(x) - (21*x*(Common::pow4(x) - 60*Common::pow2(x) + 495))*cos(x))/Common::pow7(x);
}

double Common::SphericalBesselJ8(double x) {
    if(fabs(x) < 1e-4)
        return Common::pow8(x)/34459425 - Common::pow10(x)/1309458150;
    else
        return ((Common::pow8(x) - 630*Common::pow6(x) + 51975*Common::pow4(x) - 945945*Common::pow2(x) + 2027025)*sin(x) + (9*x*(4*Common::pow6(x) - 770*Common::pow4(x) + 30030*Common::pow2(x) - 225225))*cos(x))/Common::pow9(x);
}


void Common::write(FILE* stream, const char* format, ...) {
    va_list args;
    char dest[6000];
    va_start(args, format);
    vsnprintf(dest, 2048, format, args);
    va_end(args);
    fprintf(stream, dest);
}


void Common::info(const char* format, ...) {
    va_list ap;
    va_start(ap, format);
    vfprintf(stdout, format, ap);
    va_end(ap);
    fflush(stdout);
}

void Common::verbose(const char* format, ...) {
    va_list ap;
    va_start(ap, format);
#ifdef VERBOSE
    vfprintf(stdout, format, ap);
    va_end(ap);
    fflush(stdout);
#endif
}

void Common::debug(const char* format, ...) {
    va_list ap;
    va_start(ap, format);
#ifdef DEBUG
    vfprintf(stdout, format, ap);
    va_end(ap);
    fflush(stdout);
#endif
}

void Common::warning(const char* format, ...) {
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    fflush(stderr);
}

void Common::error(const char* format, ...) {
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    fflush(stderr);
    abort();
}


void Common::throw_error(const char* msg, std::string file, int lineno) 
{
    std::string emsg(msg);
    std::string out = emsg + ", file: " __FILE__ + ", line: " + std::to_string(__LINE__);
    throw std::runtime_error(out);  
}

void Common::throw_error(std::string msg, std::string file, int lineno) 
{
    std::string out = msg + " (file: " __FILE__ + ", line: " + std::to_string(__LINE__) + ")";
    throw std::runtime_error(out);  
}
