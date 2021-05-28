#pragma once

class Communication {
  public:

    Communication(/*neighbour */);
    
   /****************
   * Initialize communication:
   * ******************/
    static void init_parallel(/*args = num of processors?*/);

    /****************
    * Finalize communication:
    * ******************/
    static void finalize();

    /****************
    * Communicate a field:
    * ******************/
    static void communicate(/* field = field class? or vectors of subfields? */);

    /*******************
     * Find a minimum value across all processes: //adaptive time step
     * *******************/
    static double reduce_min(/* value */ );

      /*******************
     * Find a minimum value across all processes: //residual
     * *******************/
    static double reduce_sum(/*value */ );

};