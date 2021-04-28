#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {

const std::vector<double> u_fixedWallBoundary_values(_u.jmax,0); 
const std::vector<double> v_fixedWallBoundary_values(_u.jmax,0); 

_U.set_col(u_fixedWallBoundary_values, jmax+1);




}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {

const std::vector<double>u_movingwallboundry_values (u.imax(),1);// fix U vs u
const std::vector<double>v_movingwallboundry_values (_v.imax(),0);
    _u.set_row(u_movingwallboundry_values, imax+1);
    _v.set_row(u_movingwallboundry_values, imax+1);
}
