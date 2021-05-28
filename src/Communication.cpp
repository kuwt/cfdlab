#include "Communication.hpp"

Communication::Communication(/*neighbour */)
{

    
}

/****************
 * Initialize communication:
 * ******************/
void Communication::init_parallel(/*args = num of processors?*/)
{

}

/****************
* Finalize communication:
* ******************/
void Communication::finalize()
{

}

/****************
* Communicate a field:
* ******************/
void Communication::communicate(/* field = field class? or vectors of subfields? */)
{

}

/*******************
 * Find a minimum value across all processes: //adaptive time step
 * *******************/
double Communication::reduce_min(/* value */ )
{

}

    /*******************
 * Find a minimum value across all processes: //residual
 * *******************/
double Communication::reduce_sum(/*value */ )
{

}

