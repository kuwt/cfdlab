#include "Fields.hpp"

#include <algorithm>
#include <iostream>

// Temporarily no need to do anything for the initial velocity unless there is a non zero initial velocity config
Fields::Fields(const Grid &grid, double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI,
               double GX, double GY, bool energy_on = false, double TI = 0.0, double Pr = 1.0, double beta = 0.1)
    : _nu(nu), _dt(dt), _tau(tau), _gx(GX), _gy(GY), _energy_on(energy_on), _Pr(Pr), _beta(beta) {
    _U = Matrix<double>(imax + 2, jmax + 2, 0.0); // NOTE: construct for the whole domain(including boundary)
    _V = Matrix<double>(imax + 2, jmax + 2, 0.0); // NOTE: UI, VI, PI: initial condition
    _P = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _T = Matrix<double>(imax + 2, jmax + 2, 0.0);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);

    // init the fields only at fluid cells
    for (auto fluid_cell : grid.fluid_cells()) {
        int i = fluid_cell->i();
        int j = fluid_cell->j();
        _U(i, j) = UI;
        _V(i, j) = VI;
        _P(i, j) = PI;
        _T(i, j) = TI;
    }
}

void Fields::calculate_temperature(Grid &grid) {
    if (_energy_on) {
        Matrix<double> _TNext = _T;
        double diffusion_term, convective_term;
        double alpha = _nu / _Pr; // thermal diffusivity
        for (auto fluid_cell : grid.fluid_cells()) {
             if (fluid_cell->isGhost()){
                continue;
            }
            int i = fluid_cell->i();
            int j = fluid_cell->j();

            convective_term = Discretization::convection_T(_T, _U, _V, i, j);
            diffusion_term = Discretization::diffusion(_T, i, j);

            _TNext(i, j) = _T(i, j) + _dt * (alpha * diffusion_term - convective_term);
        }
        _T = _TNext;
    }
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
    for (auto fluid_cell : grid.fluid_cells()) {
        if (fluid_cell->isGhost()){
            continue;
        }
        int i = fluid_cell->i();
        int j = fluid_cell->j();

        convective_term = Discretization::convection_u(_U, _V, i, j);
        diffusion_term = Discretization::diffusion(_U, i, j);
        _F(i, j) = _U(i, j) + _dt * (_nu * diffusion_term - convective_term);
        if (_energy_on) {
            _F(i, j) +=  -0.5 * _beta * _dt * _gx * (_T(i, j) + _T(i + 1, j));
        } else {
            _F(i, j) +=  _dt * _gx;
        }

        convective_term = Discretization::convection_v(_U, _V, i, j);
        diffusion_term = Discretization::diffusion(_V, i, j);
         _G(i, j) = _V(i, j) + _dt * (_nu * diffusion_term - convective_term);
        if (_energy_on) {
            _G(i, j) += -0.5 * _beta * _dt * _gy * (_T(i, j) + _T(i, j + 1));
        } else {
            _G(i, j) += _dt * _gy;
        }
    }

    // calculate for boundary cells
    for (auto boundary_cell : grid.fixed_wall_cells()) {
        if (boundary_cell->isGhost()){
            continue;
        }
        int i = boundary_cell->i();
        int j = boundary_cell->j();

        if (boundary_cell->is_border(border_position::TOP) &&
            boundary_cell->is_border(border_position::RIGHT)) // B_NE cells
        {
            _F(i, j) = _U(i, j);
            _G(i, j) = _V(i, j);
        } else if (boundary_cell->is_border(border_position::TOP) &&
                   boundary_cell->is_border(border_position::LEFT)) // B_NW cells
        {
            _F(i - 1, j) = _U(i - 1, j);
            _G(i, j) = _V(i, j);
        } else if (boundary_cell->is_border(border_position::BOTTOM) &&
                   boundary_cell->is_border(border_position::RIGHT)) // B_SE cells
        {
            _F(i, j) = _U(i, j);
            _G(i, j - 1) = _V(i, j - 1);
        } else if (boundary_cell->is_border(border_position::BOTTOM) &&
                   boundary_cell->is_border(border_position::LEFT)) // B_SW cells
        {
            _F(i - 1, j) = _U(i - 1, j);
            _G(i, j - 1) = _V(i, j - 1);
        } else if (boundary_cell->is_border(border_position::RIGHT)) // B_E cells
        {
            _F(i, j) = _U(i, j);
        } else if (boundary_cell->is_border(border_position::LEFT)) // B_W cells
        {
            _F(i - 1, j) = _U(i - 1, j);
        } else if (boundary_cell->is_border(border_position::TOP)) // B_N cells
        {
            _G(i, j) = _V(i, j);
        } else if (boundary_cell->is_border(border_position::BOTTOM)) // B_S cells
        {
            _G(i, j - 1) = _V(i, j - 1);
        }
    }

    for (auto inflow_cell : grid.inflow_cells()) {
        if (inflow_cell->isGhost()){
            continue;
        }
        // assume inflow from left
        int i = inflow_cell->i();
        int j = inflow_cell->j();
        _F(i, j) = _U(i, j);
    }
    for (auto outflow_cell : grid.outflow_cells()) {
        if (outflow_cell->isGhost()){
            continue;
        }
        // assume outflow to right
        int i = outflow_cell->i();
        int j = outflow_cell->j();
        _F(i - 1, j) = _U(i - 1, j);
    }
    for (auto moving_cell : grid.moving_wall_cells()) {
        if (moving_cell->isGhost()){
            continue;
        }
        // assume moving wall only on top
        int i = moving_cell->i();
        int j = moving_cell->j();
        _G(i, j - 1) = _V(i, j - 1);
    }
}

void Fields::calculate_rs(Grid &grid) {
    // Equation(11)
    int imax = grid.imax();
    int jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy();

    for (auto fluid_cell : grid.fluid_cells()) {
        if (fluid_cell->isGhost()){
            continue;
        }
        int i = fluid_cell->i();
        int j = fluid_cell->j();
        _RS(i, j) = 1 / _dt * ((_F(i, j) - _F(i - 1, j)) / dx + (_G(i, j) - _G(i, j - 1)) / dy);
    }
}

void Fields::calculate_velocities(Grid &grid) {
    // Equation(7)(8) after solving p at time n+1
    int imax = grid.imax();
    int jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy();

    for (auto fluid_cell : grid.fluid_cells()) {
        if (fluid_cell->isGhost()){
            continue;
        }
        int i = fluid_cell->i();
        int j = fluid_cell->j();
        if (fluid_cell->neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
            _U(i, j) = _F(i, j) - _dt / dx * (_P(i + 1, j) - _P(i, j));
        }

        if (fluid_cell->neighbour(border_position::TOP)->type() == cell_type::FLUID) {
            _V(i, j) = _G(i, j) - _dt / dy * (_P(i, j + 1) - _P(i, j));
        }
    }
}

double Fields::calculate_dt(Grid &grid) {
    // stability condition Equation(13)

    if (_tau <= 0) {
        return _dt;
    }

    double umax = 0.0, vmax = 0.0;
    double dx = grid.dx();
    double dy = grid.dy();
    int imax = grid.imax();
    int jmax = grid.jmax();

    // finding umax in the domain
    // finding vmax in the domain
    for (auto fluid_cell : grid.fluid_cells()) {
        if (fluid_cell->isGhost()){
            continue;
        }
        int i = fluid_cell->i();
        int j = fluid_cell->j();
        umax = abs(_U(i, j)) > umax ? abs(_U(i, j)) : umax;
        vmax = abs(_V(i, j)) > vmax ? abs(_V(i, j)) : vmax;
    }

    // min is taken since we want dt to be smaller than all three conditions. Only the minimum will satisfy all
    // critieria.
    _dt = _tau * std::min({dx / umax, dy / vmax, 1 / (1 / (dx * dx) + 1 / (dy * dy)) / (2 * _nu)});
    if (_energy_on) {
        double alpha = _nu / _Pr;
        _dt = _tau * std::min({_dt, _tau * 1 / (1 / (dx * dx) + 1 / (dy * dy)) / (2 * alpha)});
   }
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

bool Fields::energy_on() { return _energy_on; }


Matrix<double> &Fields::u_matrix() { return _U; }
Matrix<double> &Fields::v_matrix() { return _V; }
Matrix<double> &Fields::F_matrix() { return _F; }
Matrix<double> &Fields::G_matrix() { return _G; }
Matrix<double> &Fields::T_matrix() { return _T; }
