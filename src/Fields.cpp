#include "Fields.hpp"

#include <algorithm>
#include <iostream>


Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI)
    : _nu(nu), _dt(dt), _tau(tau) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI); // NOTE: construct for the whole domain(including boundary)
    _V = Matrix<double>(imax + 2, jmax + 2, VI); // NOTE: UI, VI, PI: initial condition
    _P = Matrix<double>(imax + 2, jmax + 2, PI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

void Fields::calculate_fluxes(Grid &grid) {
    /* NOTE:
     inside grid, there are different cell types: fluid, fixed_wall, moving_wall
     for fluid cells(inner cells) Equation (9)(10)
     for boundary cells Equation (17)
     In fields, we don't know the size of domain(but can calculate from 
     size of matrix _U actually)
    */
    
    int imax = grid.imax();
    int jmax = grid.jmax();
    double diffusion_term, convective_term;
    //[0,1,...imax, imax+1][0,1,...,jmax,jmax+1]
    
    // calculate for fluid_cells(inner cells)
    // NOTE: ranges of i, j for _F,_G are different
    for (int i = 1; i <= imax-1; ++i){ 
        for (int j = 1; j <= jmax; ++j){
            convective_term = Discretization::convection_u(_U, _V, i, j);
            diffusion_term = Discretization::diffusion(_U, i, j);           
            _F(i, j) = _U(i, j) + _dt * (_nu * diffusion_term - convective_term + _gx);           
        }
    }

    for (int i = 1; i <= imax; ++i){ 
        for (int j = 1; j <= jmax-1; ++j){
            convective_term = Discretization::convection_v(_U, _V, i, j);
            diffusion_term = Discretization::diffusion(_V, i, j);            
            _G(i, j) = _V(i, j) + _dt * (_nu * diffusion_term - convective_term + _gy);
        }
    }

    // calculate for boundary cells
    for (int j = 1; j <= jmax; ++j){
        _F(0, j) = _U(0, j);
        _F(imax, j) = _U(imax, j);
    } 

    for (int i = 1; i <= imax; ++i){
        _G(i, 0) = _V(i, 0);
        _G(i, jmax) = _V(i, jmax);
    }


} //TODO: //DONE

void Fields::calculate_rs(Grid &grid) {
    // Equation(11)
    int imax = grid.imax();
    int jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy();
    for (int i = 1; i <= imax; ++i){
        for (int j = 1; j <= jmax; ++j){
            _RS(i,j) = 1 / _dt * ((_F(i,j) - _F(i-1, j))/dx + \
                                  (_G(i,j) - _G(i, j-1))/dy);
        }
    }
} //TODO: //DONE

void Fields::calculate_velocities(Grid &grid) {
    // Equation(7)(8) after solving p at time n+1
    int imax = grid.imax();
    int jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy();

    for (int i = 1; i <= imax-1; ++i){ 
        for (int j = 1; j <= jmax; ++j){
           _U(i, j) = _F(i, j) - _dt/dx * (_P(i+1, j) - _P(i, j));           
        }
    }

    for (int i = 1; i <= imax; ++i){ 
        for (int j = 1; j <= jmax-1; ++j){
            _V(i, j) = _G(i, j) - _dt/dy * (_P(i, j+1) - _P(i, j)); 
        }
    }
} //TODO: //DONE

double Fields::calculate_dt(Grid &grid) { 
    // stability condition Equation(13)

    if (_tau <= 0){
        return _dt; 
    } 
    
    double umax = 0.0, vmax = 0.0;
    double dx = grid.dx();
    double dy = grid.dy();
    int imax = grid.imax();
    int jmax = grid.jmax();

    for (int i = 0; i <= imax; ++i){
        for (int j = 0; j<= jmax+1; ++j){
            umax =  _U(i, j) > umax ? _U(i,j) : umax;
        }
    }
    for (int i = 0; i <= imax+1; ++i){
        for (int j = 0; j<= jmax; ++j){
            vmax =  _V(i, j) > vmax ? _V(i,j) : vmax;
        }
    } 
    // umax = * std::max_element(_U.data(), _U.data()+_U.size());
    // vmax = * std::max_element(_V.data(), _U.data()+_V.size());
    _dt = _tau * std::min( {dx/abs(umax), 
                            dy/abs(vmax), 
                            1/(1/(dx*dx) + 1/(dy*dy)) /(2*_nu)});
    return _dt; 
} //TODO i,j range not sure

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; } 

double Fields::dt() const { return _dt; }