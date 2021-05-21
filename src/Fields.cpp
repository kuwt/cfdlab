#include "Fields.hpp"

#include <algorithm>
#include <iostream>

// Temporarily no need to do anything for the initial velocity unless there is a non zero initial velocity config
Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI, bool energy_on, double TI=0.0)
    : _nu(nu), _dt(dt), _tau(tau) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI); // NOTE: construct for the whole domain(including boundary)
    _V = Matrix<double>(imax + 2, jmax + 2, VI); // NOTE: UI, VI, PI: initial condition
    _P = Matrix<double>(imax + 2, jmax + 2, PI);
    _T = Matrix<double>(imax + 2, jmax + 2, TI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _energy_on = energy_on;
}

// need to use the grid information more extensively by using for each cell, checking if the cell is a fluid or not
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
    // implementation of ws1 equation (9)
    // implementation of ws1 equation (10)
    for (auto fluid_cell : grid.fluid_cells()){
        int i = fluid_cell->i();
        int j = fluid_cell->j();

        convective_term = Discretization::convection_u(_U, _V, i, j);
        diffusion_term = Discretization::diffusion(_U, i, j);           
        _F(i, j) = _U(i, j) + _dt * (_nu * diffusion_term - convective_term + _gx);
    
        convective_term = Discretization::convection_v(_U, _V, i, j);
        diffusion_term = Discretization::diffusion(_V, i, j);            
        _G(i, j) = _V(i, j) + _dt * (_nu * diffusion_term - convective_term + _gy);
    }


    // calculate for boundary cells
    for (auto boundary_cell : grid.fixed_wall_cells()){
        int i = boundary_cell->i();
        int j = boundary_cell->j();

        if (boundary_cell->is_border(border_position::TOP) && boundary_cell->is_border(border_position::RIGHT) )  // B_NE cells
        {
            _F(i, j) = _U(i, j);
            _G(i, j) = _V(i, j);
        }
        else if (boundary_cell->is_border(border_position::TOP) && boundary_cell->is_border(border_position::LEFT) )  // B_NW cells
        {
            _F(i-1, j) = _U(i-1, j);
            _G(i, j) = _V(i, j);
        }
        else if (boundary_cell->is_border(border_position::BOTTOM) && boundary_cell->is_border(border_position::RIGHT) )  // B_SE cells
        {
            _F(i, j) = _U(i, j);
            _G(i, j-1) = _V(i, j-1);
        }
        else if (boundary_cell->is_border(border_position::BOTTOM) && boundary_cell->is_border(border_position::LEFT) )  // B_SW cells
        {
            _F(i-1, j) = _U(i-1, j);
            _G(i, j-1) = _V(i, j-1);
        }   
        else if (boundary_cell->is_border(border_position::RIGHT) ) // B_E cells
        {
            _F(i, j) = _U(i, j);
        }
        else if (boundary_cell->is_border(border_position::LEFT) ) // B_W cells
        {
            _F(i-1, j) = _U(i-1, j);
        }
        else if (boundary_cell->is_border(border_position::TOP) )  // B_N cells
        {
            _G(i, j) = _V(i, j);
        }
        else if (boundary_cell->is_border(border_position::BOTTOM) )  // B_S cells
        {
            _G(i, j-1) = _V(i, j-1);
        }
        
    }

    for (auto inflow_cell : grid.inflow_cells()){
        // assume inflow from left
        int i = inflow_cell->i();
        int j = inflow_cell->j();
        _F(i , j) = _U(i, j);

    }
    for (auto outflow_cell : grid.outflow_cells()){
        // assume outflow to right
        int i = outflow_cell->i();
        int j = outflow_cell->j();
        _F(i-1 , j) = _U(i-1, j);
    }
    for (auto moving_cell : grid.moving_wall_cells()){
        // assume moving wall only on top
        int i = moving_cell->i();
        int j = moving_cell->j();
        _G(i, j-1) = _V(i, j-1);
    }

}

void Fields::calculate_rs(Grid &grid) {
    // Equation(11)
    int imax = grid.imax();
    int jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy();

    for (auto fluid_cell : grid.fluid_cells()){
        int i = fluid_cell->i();
        int j = fluid_cell->j();
        _RS(i,j) = 1 / _dt * ((_F(i,j) - _F(i-1, j))/dx + \
                                  (_G(i,j) - _G(i, j-1))/dy);
    }
}

void Fields::calculate_velocities(Grid &grid) {
    // Equation(7)(8) after solving p at time n+1
    int imax = grid.imax();
    int jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy();

    for (auto fluid_cell : grid.fluid_cells()){
        int i = fluid_cell->i();
        int j = fluid_cell->j();
        if (fluid_cell->neighbour(border_position::RIGHT)->type() == cell_type::FLUID)
        {
            _U(i, j) = _F(i, j) - _dt/dx * (_P(i+1, j) - _P(i, j));    
        }
      
        if (fluid_cell->neighbour(border_position::TOP)->type() == cell_type::FLUID)
        {
            _V(i, j) = _G(i, j) - _dt/dy * (_P(i, j+1) - _P(i, j)); 
        }
    }
} 

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

    // finding umax in the domain
    for (int i = 0; i <= imax; ++i){
        for (int j = 0; j<= jmax; ++j){
            umax =  _U(i, j) > umax ? _U(i,j) : umax;
        }
    }

    // finding vmax in the domain
    for (int i = 0; i <= imax; ++i){
        for (int j = 0; j<= jmax; ++j){
            vmax =  _V(i, j) > vmax ? _V(i,j) : vmax;
        }
    } 
    // umax = * std::max_element(_U.data(), _U.data()+_U.size());
    // vmax = * std::max_element(_V.data(), _U.data()+_V.size());

    // min is taken since we want dt to be smaller than all three conditions. Only the minimum will satisfy all critieria.
    _dt = _tau * std::min( {dx/abs(umax), 
                            dy/abs(vmax), 
                            1/(1/(dx*dx) + 1/(dy*dy)) /(2*_nu),
                            1/(1/(dx*dx) + 1/(dy*dy)) /(2*_Pr)});

    return _dt; 
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; } 

double Fields::dt() const { return _dt; }

double &Fields::T(int i, int j) { return _T(i, j); }