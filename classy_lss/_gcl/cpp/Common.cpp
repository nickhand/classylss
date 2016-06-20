#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <stdexcept>

#include "Common.h"

double Common::SphericalBesselJ0(double x) {
    if(fabs(x) < 1e-4)
        return 1 - pow2(x)/6 + pow4(x)/120 - pow6(x)/5040;
    else
        return sin(x)/x;
}

double Common::SphericalBesselJ1(double x) {
    if(fabs(x) < 1e-4)
        return x/3 - pow3(x)/30 + pow5(x)/840 - pow7(x)/45360;
    else
        return (sin(x) - x*cos(x))/pow2(x);
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
