#include "Boundary.hpp"
#include <cmath>
#include <iostream>

static const double zeroEpilon = 1e-05;

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {
    for (int cell_iter = 0; cell_iter < _cells.size(); ++cell_iter) {
        auto pcell = _cells[cell_iter];
        int i = pcell->i();
        int j = pcell->j();

        double walltemperature = 0;
        if (_wall_temperature.find(pcell->wall_id()) != _wall_temperature.end()) {
            walltemperature = _wall_temperature[pcell->wall_id()];
        }

        if (pcell->is_border(border_position::TOP) && pcell->is_border(border_position::RIGHT)) // B_NE cells
        {
            field.u(i, j) = 0;
            field.v(i, j) = 0;
            field.u(i - 1, j) = -field.u(i - 1, j + 1);
            field.v(i - 1, j) = -field.v(i + 1, j - 1);

            if (walltemperature < (-1 + zeroEpilon)) { // Neumann adiabatic
                field.T(i, j) = (field.T(i + 1, j) + field.T(i, j + 1)) / 2;
            } else { // Dirichlet
                field.T(i, j) = 2 * walltemperature - (field.T(i + 1, j) + field.T(i, j + 1)) / 2;
            }
        }

        else if (pcell->is_border(border_position::TOP) && pcell->is_border(border_position::LEFT)) // B_NW cells
        {
            field.v(i, j) = 0;
            field.u(i - 1, j) = 0;
            field.u(i, j) = -field.u(i, j + 1);
            field.v(i, j - 1) = -field.v(i - 1, j - 1);

            if (walltemperature < (-1 + zeroEpilon)) { // Neumann adiabatic
                field.T(i, j) = (field.T(i - 1, j) + field.T(i, j + 1)) / 2;
            } else { // Dirichlet
                field.T(i, j) = 2 * walltemperature - (field.T(i - 1, j) + field.T(i, j + 1)) / 2;
            }

        } else if (pcell->is_border(border_position::BOTTOM) && pcell->is_border(border_position::RIGHT)) // B_SE cells
        {
            field.u(i, j) = 0;
            field.v(i, j - 1) = 0;
            field.u(i - 1, j) = -field.u(i - 1, j - 1);
            field.v(i, j) = -field.v(i + 1, j);

            if (walltemperature < (-1 + zeroEpilon)) { // Neumann adiabatic
                field.T(i, j) = (field.T(i + 1, j) + field.T(i, j - 1)) / 2;
            } else { // Dirichlet
                field.T(i, j) = 2 * walltemperature - (field.T(i + 1, j) + field.T(i, j - 1)) / 2;
            }
        } else if (pcell->is_border(border_position::BOTTOM) && pcell->is_border(border_position::LEFT)) // B_SW cells
        {
            field.u(i - 1, j) = 0;
            field.v(i, j - 1) = 0;
            field.u(i, j) = -field.u(i, j - 1);
            field.v(i, j) = -field.v(i - 1, j);

            if (walltemperature < (-1 + zeroEpilon)) { // Neumann adiabatic
                field.T(i, j) = (field.T(i - 1, j) + field.T(i, j - 1)) / 2;
            } else { // Dirichlet
                field.T(i, j) = 2 * walltemperature - (field.T(i - 1, j) + field.T(i, j - 1)) / 2;
            }
        } else if (pcell->is_border(border_position::RIGHT)) // B_E cells
        {
            field.u(i, j) = 0;
            field.v(i, j) = -field.v(i + 1, j);
            field.v(i, j - 1) = -field.v(i + 1, j - 1);

            if (walltemperature < (-1 + zeroEpilon)) { // Neumann adiabatic
                field.T(i, j) = field.T(i + 1, j);
            } else { // Dirichlet
                field.T(i, j) = 2 * walltemperature - field.T(i + 1, j);
            }
        } else if (pcell->is_border(border_position::LEFT)) // B_W cells
        {
            field.u(i - 1, j) = 0;
            field.v(i, j - 1) = -field.v(i - 1, j - 1);
            field.v(i, j) = -field.v(i - 1, j);

            if (walltemperature < (-1 + zeroEpilon)) { // Neumann adiabatic
                field.T(i, j) = field.T(i - 1, j);
            } else { // Dirichlet
                field.T(i, j) = 2 * walltemperature - field.T(i - 1, j);
            }
        } else if (pcell->is_border(border_position::TOP)) // B_N cells
        {
            field.v(i, j) = 0;
            field.u(i - 1, j) = -field.u(i - 1, j + 1);
            field.u(i, j) = -field.u(i, j + 1);

            if (walltemperature < (-1 + zeroEpilon)) { // Neumann adiabatic
                field.T(i, j) = field.T(i, j + 1);
            } else { // Dirichlet
                field.T(i, j) = 2 * walltemperature - field.T(i, j + 1);
            }

        } else if (pcell->is_border(border_position::BOTTOM)) // B_S cells
        {
            field.v(i, j - 1) = 0;
            field.u(i, j) = -field.u(i, j - 1);
            field.u(i - 1, j) = -field.u(i - 1, j - 1);

            if (walltemperature < (-1 + zeroEpilon)) { // Neumann adiabatic
                field.T(i, j) = field.T(i, j - 1);
            } else { // Dirichlet
                field.T(i, j) = 2 * walltemperature - field.T(i, j - 1);
            }
        }
    }
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {
    /*********
     *  Important warning! Currently only works for Lid Driven Cavity Case
     * *******/
    for (int cell_iter = 0; cell_iter < _cells.size(); ++cell_iter) {
        auto pcell = _cells[cell_iter];
        int i = pcell->i();
        int j = pcell->j();

        if (pcell->is_border(border_position::BOTTOM)) // top wall
        {
            field.v(i, j - 1) = 0; // the minus 1 is to account the fact the the top border is just the top ghost border
            field.u(i, j) = 2 * _wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(i, j - 1);
        }
    }
}

InFlowBoundary::InFlowBoundary(std::vector<Cell *> cells, double inflow_velocity_x, double inflow_velocity_y)
    : _cells(cells), _inflow_velocity_x(inflow_velocity_x), _inflow_velocity_y(inflow_velocity_y) {}

void InFlowBoundary::apply(Fields &field) {
    /*********
     *  Important warning! Currently only works When inflow is at the left
     *  No temperature effect is considered here too.
     * *******/

    // Dirichlet BC for inflow
    for (int cell_iter = 0; cell_iter < _cells.size(); ++cell_iter) {
        auto pcell = _cells[cell_iter];
        int i = pcell->i();
        int j = pcell->j();

        // assume inflow from left
        field.u(i, j) = _inflow_velocity_x;
        field.v(i, j) = 2 * _inflow_velocity_y - field.v(i + 1, j);
        field.v(i, j - 1) = 2 * _inflow_velocity_y - field.v(i + 1, j - 1);
    }
}

OutFlowBoundary::OutFlowBoundary(std::vector<Cell *> cells) : _cells(cells) {}

void OutFlowBoundary::apply(Fields &field) {
    /*********
     *  Important warning! Currently only works When outflow is at the right
     * No temperature effect is considered here too.
     * *******/

    // Neumann BC for outflow(free-outflow)
    for (int cell_iter = 0; cell_iter < _cells.size(); ++cell_iter) {
        auto pcell = _cells[cell_iter];
        int i = pcell->i();
        int j = pcell->j();

        // assume outflow to right
        // i() - 1 is the exact position of the boundary and we need to make it equal to i() -2
        // i() never appears in the simulation
        field.u(i - 1, j) = field.u(i - 2, j);
        field.v(i - 1, j) = field.v(i - 2, j);
        field.v(i - 1, j - 1) = field.v(i - 2, j - 1);
    }
}
