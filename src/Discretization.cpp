#include "Discretization.hpp"

#include <cmath>

double Discretization::_dx = 0.0;
double Discretization::_dy = 0.0;
double Discretization::_gamma = 0.0;

Discretization::Discretization(double dx, double dy, double gamma) {
    _dx = dx;
    _dy = dy;
    _gamma = gamma;
}

double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    // du2/dx + duv/dy
    // implementation of equation (4)
    double du2dx = 1/_dx * ( pow((U(i,j)+U(i+1,j))/2, 2) - \
                             pow((U(i-1,j)+U(i,j))/2, 2) )+ \
                   _gamma / _dx * ( fabs(U(i,j)+U(i+1,j)) * (U(i,j)-U(i+1,j)) / 4 - \
                                    fabs(U(i-1,j)+U(i,j)) * (U(i-1,j)-U(i,j)) / 4);
    double duvdy = 1/_dy * ( (V(i,j)+V(i+1,j)) * (U(i,j)+U(i,j+1)) / 4 - \
                             (V(i,j-1)+V(i+1,j-1)) * (U(i,j-1)+U(i,j)) / 4) + \
                   _gamma / _dy * ( fabs(V(i,j)+V(i+1,j)) * (U(i,j)-U(i,j+1)) / 4 - \
                                    fabs(V(i,j-1)+V(i+1,j-1)) * (U(i,j-1)-U(i,j)) / 4);
    return du2dx + duvdy;
                             
}//TODO: //DONE

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    // duv/_dx + dv2/dy
    double duvdx = 1/_dx * ( (U(i,j)+U(i,j+1)) * (V(i,j)+V(i+1,j)) / 4 - \
                             (U(i-1,j)+U(i-1,j+1)) * (V(i-1,j)+V(i,j)) / 4) + \
                   _gamma / _dx * ( fabs(U(i,j)+U(i,j+1)) * (V(i,j)-V(i+1,j)) / 4 - \
                                    fabs(U(i-1,j)+U(i-1,j+1)) * (V(i-1,j)-V(i,j)) / 4);
    double dv2dy = 1/_dy * ( pow((V(i,j) + V(i,j+1))/2, 2) - \
                             pow((V(i,j-1) + V(i,j))/2, 2) )+ \
                   _gamma / _dy * ( fabs(V(i,j)+V(i,j+1)) * (V(i,j)-V(i,j+1)) / 4 - \
                                    fabs(V(i,j-1)+V(i,j)) * (V(i,j-1)-V(i,j)) / 4);
    return duvdx + dv2dy;

    
}//TODO: //DONE

double Discretization::diffusion(const Matrix<double> &A, int i, int j) {
    // U: nu * (d2u/dx2 + d2u/dy2)
    // V: nu * (d2v/dy2 + d2v/dy2)
    // the same operation applied for both U and V: U ==> A
    double result = ( A(i+1,j) - 2.0*A(i,j) + A(i-1,j) ) / (_dx * _dx) +
                    ( A(i,j+1) - 2.0*A(i,j) + A(i,j-1) ) / (_dx * _dx);
    return result;
}//TODO: //DONE

double Discretization::laplacian(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (_dx * _dx) +
                    (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (_dx * _dx);
    return result;
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {}//TODO: