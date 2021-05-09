#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {} 

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {} //TODO: todo later, currently no temperature

void FixedWallBoundary::apply(Fields &field) {
    for (int cell_iter = 0; cell_iter < _cells.size(); ++cell_iter)
    {
        // set well defined point
        if (_cells[cell_iter]->is_border(border_position::RIGHT) ) // left wall means the right handside is the real border
        {
             field.u(_cells[cell_iter]->i(),_cells[cell_iter]->j()) = 0;
        }
        if (_cells[cell_iter]->is_border(border_position::LEFT) )  // right wall
        {
             field.u(_cells[cell_iter]->i()-1,_cells[cell_iter]->j()) = 0; //the minus 1 is to account the fact the the right border is just the right ghost border
        }
        if (_cells[cell_iter]->is_border(border_position::TOP) )  // bottom wall
        {
             field.v(_cells[cell_iter]->i(),_cells[cell_iter]->j()) = 0;
        }

         // set intermediate point
        if (_cells[cell_iter]->is_border(border_position::RIGHT) ) // left wall
        {
             field.v(_cells[cell_iter]->i(),_cells[cell_iter]->j()) = -field.v(_cells[cell_iter]->i()+1,_cells[cell_iter]->j());
        }
        if (_cells[cell_iter]->is_border(border_position::LEFT) )  // right wall
        {
              field.v(_cells[cell_iter]->i(),_cells[cell_iter]->j()) = -field.v(_cells[cell_iter]->i()-1,_cells[cell_iter]->j());
        }
        if (_cells[cell_iter]->is_border(border_position::TOP) )  // bottom wall
        {
             field.u(_cells[cell_iter]->i(),_cells[cell_iter]->j()) =  -field.u(_cells[cell_iter]->i(),_cells[cell_iter]->j()+1);
        }
    }
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}//TODO: todo later, currently no temperature

void MovingWallBoundary::apply(Fields &field) {
    
    for (int cell_iter = 0; cell_iter < _cells.size(); ++cell_iter)
    {
        if (_cells[cell_iter]->is_border(border_position::BOTTOM) ) // top wall
        {
            field.v(_cells[cell_iter]->i(),_cells[cell_iter]->j()-1) = 0; //the minus 1 is to account the fact the the top border is just the top ghost border
            field.u(_cells[cell_iter]->i(),_cells[cell_iter]->j()) = 2*_wall_velocity[LidDrivenCavity::moving_wall_id] -field.u(_cells[cell_iter]->i(),_cells[cell_iter]->j()-1);
        }
    }
}



InFlowBoundary::InFlowBoundary(std::vector<Cell *> cells, double inflow_velocity_x, double inflow_velocity_y)
    : _cells(cells), _inflow_velocity_x(inflow_velocity_x), _inflow_velocity_y(inflow_velocity_y) {}

void InFlowBoundary::apply(Fields &field) {
    // to do:  diriclet
}


OutFlowBoundary::OutFlowBoundary(std::vector<Cell *> cells)
    : _cells(cells){}

void OutFlowBoundary::apply(Fields &field) {
    // to do: neuman
}