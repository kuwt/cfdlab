#include "Communication.hpp"
#include <mpi.h>

Communication::Communication()
{}

/****************
 * Initialize communication:
 * ******************/
void Communication::init_parallel(int argc, char ** argv, int &rank, int &num_proc)
{
    int mpi_status = MPI_SUCCESS;
    rank = -1;
    num_proc = -1;
    mpi_status = MPI_Init(&argc, &argv);
    if (MPI_SUCCESS != mpi_status){
        std::cerr << "MPI_Init fail.\n";
    }
    mpi_status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (MPI_SUCCESS != mpi_status){
        std::cerr << "MPI_Comm_rank fail.\n";
    }
    mpi_status = MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    if (MPI_SUCCESS != mpi_status){
        std::cerr << "MPI_Comm_size fail.\n";
    }

    if (rank == 0){
        std::cout << "My rank is " << rank << std::endl;
        std::cout << "There are  " << num_proc << "processes."<< std::endl;
    }
}

/****************
* Finalize communication:
* ******************/
void Communication::finalize()
{
    int mpi_status = MPI_SUCCESS;
    int rank = -1;
    mpi_status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (MPI_SUCCESS != mpi_status){
        std::cerr << "MPI_Comm_rank fail.\n";
    }

    mpi_status = MPI_Barrier(MPI_COMM_WORLD);
    if (MPI_SUCCESS != mpi_status){
        std::cerr << "MPI_Barrier fail.\n";
    }

    mpi_status = MPI_Finalize();
    if (MPI_SUCCESS != mpi_status){
        std::cerr << "MPI_Finalize fail.\n";
    }

    if (rank == 0) {
        std::cout << "My rank is " << rank << std::endl;
        std::cout << "MPI finalized. " << std::endl;
    }
}

/****************
* Communicate a field:
* ******************/
void Communication::communicate(const Grid &grid, Matrix<double> &mat)
{
    //fill send buffer from mat
    // send to LEFT neighbor and receive from RIGHT neighbor
    // Get receive buffer to mat

    // fill send buffer from mat
    // send to RIGHT neighbor and receive from LEFT neighbor
    // Get receive buffer to mat

    // fill send buffer from mat
    // send to TOP neighbor and receive from BOTTOM neighbor
    // Get receive buffer to mat

    // fill send buffer from mat
    // send to BOTTOM neighbor and receive from TOP neighbor
    // Get receive buffer to mat
}

/*******************
 * Find a minimum value across all processes: //e.g.adaptive time step
 * *******************/
double Communication::reduce_min(const double value)
{
    double localValue = value;
    double retValue = 0;
    int mpi_status = MPI_Reduce(&localValue, &retValue, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if (MPI_SUCCESS != mpi_status){
        std::cerr << "MPI_Reduce fail.\n";
    }
    return retValue;
}

    /*******************
 * Find a minimum value across all processes: //e.g. residual
 * *******************/
double Communication::reduce_sum(const double value)
{
    double localValue = value;
    double retValue = 0;
    int mpi_status = MPI_Reduce(&localValue, &retValue, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (MPI_SUCCESS != mpi_status){
        std::cerr << "MPI_Reduce fail.\n";
    }
    return retValue;
}

