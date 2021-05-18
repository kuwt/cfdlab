#include "PressureSolver.hpp"

#include <cmath>
#include <iostream>

SOR::SOR(double omega) : _omega(omega) {}

double SOR::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {

    double dx = grid.dx();
    double dy = grid.dy();

    double coeff = _omega / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy))); // = _omega * h^2 / 4.0, if dx == dy == h

    // apply pressure boundary
 {   
        std::vector<Cell*> boundaryCells;
        boundaryCells.insert(boundaryCells.end(), grid.fixed_wall_cells().begin(),grid.fixed_wall_cells().end());
        boundaryCells.insert(boundaryCells.end(), grid.inflow_cells().begin(), grid.inflow_cells().end() );
        boundaryCells.insert(boundaryCells.end(), grid.outflow_cells().begin(), grid.outflow_cells().end() );
        boundaryCells.insert(boundaryCells.end(), grid.moving_wall_cells().begin(), grid.moving_wall_cells().end() );
        for (auto boundaryCell : boundaryCells)
        {   
            // for inflow boundary, we assume that it locates at left
            // for outflow boundary, we assume that it locates at right
            // for moving wall boundary, we assume that it locates at top
            int i = boundaryCell->i();
            int j = boundaryCell->j();
            if (boundaryCell->is_border(border_position::TOP) && boundaryCell->is_border(border_position::RIGHT) )  // B_NE cells
            {
                field.p(i, j) = (field.p(i+1, j) + field.p(i, j+1)) / 2;
            }
            else if (boundaryCell->is_border(border_position::TOP) && boundaryCell->is_border(border_position::LEFT) )  // B_NW cells
            {
                field.p(i, j) = (field.p(i-1, j) + field.p(i, j+1)) / 2;
            }
            else if (boundaryCell->is_border(border_position::BOTTOM) && boundaryCell->is_border(border_position::RIGHT) )  // B_SE cells
            {
                field.p(i, j) = (field.p(i+1, j) + field.p(i, j-1)) / 2;
            }
            else if (boundaryCell->is_border(border_position::BOTTOM) && boundaryCell->is_border(border_position::LEFT))  // B_SW cells
            {
                field.p(i, j) = (field.p(i-1, j) + field.p(i, j-1)) / 2;
            } 
            else if (boundaryCell->is_border(border_position::RIGHT) ) // B_E cells
            {
                field.p(i, j) = field.p(i+1, j);
            }
            else if (boundaryCell->is_border(border_position::LEFT) ) // B_W cells
            {
                field.p(i, j) = field.p(i-1, j);
            }
            else if (boundaryCell->is_border(border_position::TOP) )  // B_N cells
            {
                field.p(i, j) = field.p(i, j+1);
            }
            else if (boundaryCell->is_border(border_position::BOTTOM) )  // B_S cells
            {
                field.p(i, j) = field.p(i, j-1);
            }
            
        }
    }

      // update p according to stencil
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        field.p(i, j) = (1.0 - _omega) * field.p(i, j) +
                        coeff * (Discretization::sor_helper(field.p_matrix(), i, j) - field.rs(i, j));
    }

     // compute residual
    double res = 0.0;
    double rloc = 0.0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += (val * val);
    }
    {
        res = rloc / (grid.fluid_cells().size());
        res = std::sqrt(res);
    }

    return res;
}
