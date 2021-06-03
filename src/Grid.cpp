#include "Grid.hpp"
#include "Enums.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

Grid::Grid(std::string geom_name, Domain &domain) {

    _domain = domain;

    _cells = Matrix<Cell>(_domain.size_x + 2, _domain.size_y + 2);

    if (geom_name.compare("NONE")) {
        std::vector<std::vector<int>> geometry_data(_domain.domain_size_x + 2,
                                                    std::vector<int>(_domain.domain_size_y + 2, 0));
        parse_geometry_file(geom_name, geometry_data);
        assign_cell_types(geometry_data);
    } else {
        build_lid_driven_cavity();
    }
}

void Grid::build_lid_driven_cavity() {
    std::vector<std::vector<int>> geometry_data(_domain.domain_size_x + 2,
                                                std::vector<int>(_domain.domain_size_y + 2, 0));

    for (int i = 0; i < _domain.domain_size_x + 2; ++i) {
        for (int j = 0; j < _domain.domain_size_y + 2; ++j) {
            // Bottom, left and right walls: no-slip
            if (i == 0 || j == 0 || i == _domain.domain_size_x + 1) {
                geometry_data.at(i).at(j) = LidDrivenCavity::fixed_wall_id;
            }
            // Top wall: moving wall
            else if (j == _domain.domain_size_y + 1) {
                geometry_data.at(i).at(j) = LidDrivenCavity::moving_wall_id;
            }
        }
    }
    assign_cell_types(geometry_data);
}

// to do, need to adapt to geometry
// have to do the  cell neighbour assigment for each cell
void Grid::assign_cell_types(std::vector<std::vector<int>> &geometry_data) {

    int i = 0;
    int j = 0;

    for (int j_geom = _domain.jmin; j_geom < _domain.jmax; ++j_geom) {
        {
            i = 0;
        }
        for (int i_geom = _domain.imin; i_geom < _domain.imax; ++i_geom) {
            if (geometry_data.at(i_geom).at(j_geom) == 0) {
                // fluid cells
                _cells(i, j) = Cell(i, j, cell_type::FLUID);
                _fluid_cells.push_back(&_cells(i, j));

                // moving wall cells-lid Driven Cavity
            } else if (geometry_data.at(i_geom).at(j_geom) == LidDrivenCavity::moving_wall_id) {
                _cells(i, j) = Cell(i, j, cell_type::MOVING_WALL, geometry_data.at(i_geom).at(j_geom));
                _moving_wall_cells.push_back(&_cells(i, j));
                // Fixed wall Cells-Lid Driven Cavity
            } else if (geometry_data.at(i_geom).at(j_geom) == LidDrivenCavity::fixed_wall_id) {
                if (i == 0 or j == 0 or i == _domain.size_x + 1 or j == _domain.size_y + 1)
                    _cells(i, j) = Cell(i, j, cell_type::FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                _fixed_wall_cells.push_back(&_cells(i, j));

            } else if (geometry_data.at(i_geom).at(j_geom) == 1) {
                if (i == 0 or j == 0 or i == _domain.size_x + 1 or j == _domain.size_y + 1) {
                    // inflow
                    _cells(i, j) = Cell(i, j, cell_type::INFLOW, geometry_data.at(i_geom).at(j_geom));
                    _inflow_cells.push_back(&_cells(i, j));
                }
            } else if (geometry_data.at(i_geom).at(j_geom) == 2) {
                // outflow
                _cells(i, j) = Cell(i, j, cell_type::OUTFLOW, geometry_data.at(i_geom).at(j_geom));
                _outflow_cells.push_back(&_cells(i, j));

            } else if (geometry_data.at(i_geom).at(j_geom) >= 3 && geometry_data.at(i_geom).at(j_geom) <= 7) {
                // fixed wall- No Slip
                _cells(i, j) = Cell(i, j, cell_type::FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                _fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) >= 9 && geometry_data.at(i_geom).at(j_geom) <= 12) {
                // fixed wall- Free Slip
                _cells(i, j) = Cell(i, j, cell_type::FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                _fixed_wall_cells_free_slip.push_back(&_cells(i, j));
            }
            
            // ghost cell, useful for parallelism only
             if (i_geom == _domain.imin){
                _ghost_cells_left.push_back(&_cells(i, j));
             }
             if (i_geom == _domain.imax -1){
                _ghost_cells_right.push_back(&_cells(i, j));
             }
             if (i_geom == _domain.jmin){
                _ghost_cells_bottom.push_back(&_cells(i, j));
             }
              if (i_geom == _domain.jmax -1){
                _ghost_cells_top.push_back(&_cells(i, j));
             }                                                      
            ++i;
        }
        ++j;
    }

    // Iterate over all cells and assign fluid borders
    for (int i = 0; i <= _domain.size_x + 1; ++i) {
        for (int j = 0; j <= _domain.size_y + 1; ++j) {

            if (i >= 1) {
                _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
                if (_cells(i, j).type() != cell_type::FLUID) {
                    if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
                        _cells(i, j).add_border(border_position::LEFT);
                    }
                }
            }
            if (i <= _domain.size_x) {
                _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);

                if (_cells(i, j).type() != cell_type::FLUID) {
                    if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
                        _cells(i, j).add_border(border_position::RIGHT);
                    }
                }
            }
            if (j >= 1) {
                _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
                if (_cells(i, j).type() != cell_type::FLUID) {
                    if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
                        _cells(i, j).add_border(border_position::BOTTOM);
                    }
                }
            }
            if (j <= _domain.size_y) {
                _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);

                if (_cells(i, j).type() != cell_type::FLUID) {
                    if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
                        _cells(i, j).add_border(border_position::TOP);
                    }
                }
            }
        }
    }
}

void Grid::parse_geometry_file(std::string filedoc, std::vector<std::vector<int>> &geometry_data) {

    int numcols, numrows, depth;

    std::ifstream infile(filedoc);
    std::stringstream ss;
    std::string inputLine = "";

    // First line : version
    getline(infile, inputLine);
    if (inputLine.compare("P2") != 0) {
        std::cerr << "First line of the PGM file should be P2" << std::endl;
    }

    // Second line : comment
    getline(infile, inputLine);

    // Continue with a stringstream
    ss << infile.rdbuf();
    // Third line : size
    ss >> numrows >> numcols;
    // Fourth line : depth
    ss >> depth;

    int array[numrows][numcols];

    // Following lines : data
    for (int col = numcols - 1; col > -1; --col) {
        for (int row = 0; row < numrows; ++row) {
            ss >> geometry_data[row][col];
        }
    }

    infile.close();
}

int Grid::imax() const { return _domain.size_x; }
int Grid::jmax() const { return _domain.size_y; }

int Grid::imaxb() const { return _domain.size_x + 2; }
int Grid::jmaxb() const { return _domain.size_y + 2; }

Cell Grid::cell(int i, int j) const { return _cells(i, j); }

double Grid::dx() const { return _domain.dx; }

double Grid::dy() const { return _domain.dy; }

const Domain &Grid::domain() const { return _domain; }

const std::vector<Cell *> &Grid::fluid_cells() const { return _fluid_cells; }

const std::vector<Cell *> &Grid::fixed_wall_cells() const { return _fixed_wall_cells; }

const std::vector<Cell *> &Grid::moving_wall_cells() const { return _moving_wall_cells; }

const std::vector<Cell *> &Grid::inflow_cells() const { return _inflow_cells; }

const std::vector<Cell *> &Grid::outflow_cells() const { return _outflow_cells; }

const std::vector<Cell *> &Grid::fixed_wall_cells_free_slip() const { return _fixed_wall_cells_free_slip; }

const std::vector<Cell *> &Grid::ghost_cells_Left() const {return _ghost_cells_left;}

const std::vector<Cell *> &Grid::ghost_cells_Right() const {return _ghost_cells_right;}

const std::vector<Cell *> &Grid::ghost_cells_Bottom() const {return _ghost_cells_bottom;}

const std::vector<Cell *> &Grid::ghost_cells_Top() const {return _ghost_cells_top;}