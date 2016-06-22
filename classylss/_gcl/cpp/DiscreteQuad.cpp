#include "DiscreteQuad.h"

double basic_simps(const parray& y, int start, int stop, const parray& x) {

    int step(2), N(y.size());

    // slice the arrays
    parray y0 = y.slice(start, stop, step);
    parray y1 = y.slice(start+1, stop+1, step);
    parray y2 = y.slice(start+2, stop+2, step);

    // compute dx
    parray h = x.slice(1, N) - x.slice(0, -1);
    parray h0 = h.slice(start, stop, step);
    parray h1 = h.slice(start+1, stop+1, step);
    
    // the simps weights
    parray hsum = h0 + h1;
    parray hprod = h0 * h1; 
    parray h0divh1 = h0/h1;    
    parray result = (1./6.)*hsum * (y0*(2. - 1./h0divh1) + y1*hsum*hsum/hprod + y2*(2. - h0divh1));
    
    // sum
    double toret(0);
    for (double x : result) toret += x;

    return toret;
}

double TrapzIntegrate(const parray& x, const parray& y) {
    
    // trapz is just dx * bincenters
    double toret(0);    
    for (size_t i = 0; i < y.size()-1; i++)
        toret += (x[i+1] - x[i]) * 0.5 * (y[i+1] + y[i]);
    
    return toret;
}

double SimpsIntegrate(const parray& x, const parray& y) {
    
    int N(y.size());
    double toret(0.);
    double extra(0.);
    
    // straight simps requires even number of intervals
    if (N % 2 == 0) {
                
        // trapz on last point + simps on odd number
        double last_dx = x[-1] - x[-2];
        extra += 0.5*last_dx*(y[-1]+y[-2]);
        toret += basic_simps(y, 0, N-3, x); 
        
        // trapz on first point + simps on odd number
        double first_dx = x[1] - x[0];
        extra += 0.5*first_dx*(y[1]+y[0]);
        toret += basic_simps(y, 1, N-2, x);    
        
        // take the avg
        toret = 0.5*toret + 0.5*extra;
        
    } else
        toret = basic_simps(y, 0, N-2, x);
    
    return toret;
}


/******************************************************************************
 * DiscreteIntegrate
 ******************************************************************************/

double DiscreteIntegrate(int n, double* f, double h) {
    /* Use Simpson's rule when we have an even number of intervals (odd number of points) */
    if((n % 2) == 1) {
        double S = f[0] + f[n-1];
        for(int i = 1; i < n-1; i++)
            S += 2*(1 + (i%2))*f[i];
        return S*h/3;
    }
    /* Otherwise use Hollingsworth and Hunter's 3rd-order formula */
    else if(n == 2)
        return h/2 * (f[0] + f[1]);
    else if(n == 4)
        return h/8 * (3*f[0] + 9*f[1] + 9*f[2] + 3*f[3]);
    else if(n == 6)
        return h/24 * (9*f[0] + 28*f[1] + 23*f[2] + 23*f[3] + 28*f[4] + 9*f[5]);
    else {
        double S = ( 9*(f[0] + f[n-1]) + 28*(f[1] + f[n-2]) + 23*(f[2] + f[n-3]) )/24;
        for(int i = 3; i < n-3; i++)
            S += f[i];
        return S*h;
    }
}