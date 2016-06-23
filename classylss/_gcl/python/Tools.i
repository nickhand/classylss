%{
#include "Tools.h"
%}

%nodefaultctor IntegrationMethods;
%nodefaultdtor IntegrationMethods;
struct IntegrationMethods {
    enum Type {FFTLOG=0, SIMPS, TRAPZ}; 
};


parray pk_to_xi(int ell, const parray& k, const parray& pk, const parray& r, 
                double smoothing=0.5, IntegrationMethods::Type method=IntegrationMethods::FFTLOG);
parray xi_to_pk(int ell, const parray& r, const parray& xi, const parray& k, 
                double smoothing=0.005, IntegrationMethods::Type method=IntegrationMethods::FFTLOG);


