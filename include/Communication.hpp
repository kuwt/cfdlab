#pragma once

#include "Datastructures.hpp"
#include "Grid.hpp"
#include "Fields.hpp"

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
    * Communicate Domain information:
    * ******************/
    static void communicateDomainInfo(int my_rank, int their_rank, int &i_domain_min,int &i_domain_max, int &j_domain_min, int &j_domain_max);

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