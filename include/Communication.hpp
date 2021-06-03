#pragma once

#include "Datastructures.hpp"
#include "Grid.hpp"

class Communication {
  public:

    Communication();
    
   /****************
   * Initialize communication:
   * ******************/
    static void init_parallel(int argc, char ** argv, int &rank, int &num_proc);

    /****************
    * Finalize communication:
    * ******************/
    static void finalize();

    /****************
    * Communicate a field:
    * ******************/
    static void communicate(const Grid &grid, Matrix<double> &mat);

    /*******************
     * Find a minimum value across all processes: //adaptive time step
     * *******************/
    static double reduce_min(const double value);

      /*******************
     * Find a minimum value across all processes: //residual
     * *******************/
    static double reduce_sum(const double value);

};