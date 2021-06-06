#pragma once

#include "Datastructures.hpp"
#include "Grid.hpp"
#include "Fields.hpp"
#include <vector>

class Communication {
  private:
    static bool isBufInitialized;
    static std::vector<double> bufSendx;
    static std::vector<double> bufRecvx;
    static std::vector<double> bufSendy;
    static std::vector<double> bufRecvy;
    
  public:
    static int _num_proc;
    
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
    * Communicate Neighbour information:
    * ******************/
     static void communicateNeighbourInfo(int my_rank, int their_rank, int &left_neighbour_rank,int &right_neighbour_rank, int &bottom_neighbour_rank, int &top_neighbour_rank);

    /****************
    * Communicate a field:
    * ******************/
    static void communicate(const Grid &grid,
                           Matrix<double> &mat,
                           int left_neighbour_rank,
                           int right_neighbour_rank,
                           int bottom_neighbour_rank, 
                           int top_neighbour_rank);

    /*******************
     * Find a minimum value across all processes: //adaptive time step
     * *******************/
    static double reduce_min(int my_rank, const double value);

      /*******************
     * Find a minimum value across all processes: //residual
     * *******************/
    static double reduce_sum(int my_rank,  const double value);

};