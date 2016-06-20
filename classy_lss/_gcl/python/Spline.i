%{
#include "Spline.h"
%}

class Spline {
public:
    double Evaluate(double x);
    double EvaluateDerivative(double x);

    %extend {
        double __call__(double x) {
            return $self->Evaluate(x);
        }

        parray __call__(const parray& X) {
            int N = (int)X.size();
            parray Y(N);
            for(int i = 0; i < N; i++)
                Y[i] = $self->Evaluate(X[i]);
            return Y;
        }
    };
};

Spline LinearSpline(const parray& X, const parray& Y);
Spline ShiftedLinearSpline(const parray& X, const parray& Y, double tau = 0.2);
Spline CubicSpline(const parray& X, const parray& Y);
