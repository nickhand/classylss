#include "LinearPS.h"

using namespace Common;

LinearPS::LinearPS(const Cosmology& C_, double z_, double kmax_)
    : C(C_), z(z_), sigma8_z(C.Sigma8_z(z)), kmax(kmax_)
{
    
}

LinearPS::~LinearPS() 
{
    C.Clean();
}


double LinearPS::Evaluate(double k) const {
    
    double A = pow2(C.delta_H());
    double norm = pow2(sigma8_z / C.sigma8());
    return A * norm * pow(k, C.n_s())*(pow2(C.EvaluateTransfer(k)));
}
