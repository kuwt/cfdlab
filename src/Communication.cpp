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

void Communication::communicateDomainInfo(  int my_rank, 
                                            int their_rank, 
                                            int &i_domain_min,
                                            int &i_domain_max, 
                                            int &j_domain_min, 
                                            int &j_domain_max)

{
    int mpi_status = MPI_SUCCESS;
    if (my_rank == 0){   
        int sendData[4] = {i_domain_min, i_domain_max, j_domain_min,j_domain_max};
        mpi_status = MPI_Send(&sendData, 4, MPI_INT, their_rank, 0, MPI_COMM_WORLD);
        if (MPI_SUCCESS != mpi_status){
            std::cerr << "MPI_Send fail.\n";
        }
    }
    else{
        int recvData[4] = {0,0,0,0};
        MPI_Status status;
        mpi_status = MPI_Recv(&recvData, 4, MPI_INT, their_rank, 0, MPI_COMM_WORLD, &status);
        if (MPI_SUCCESS != mpi_status){
            std::cerr << "MPI_Recv fail.\n";
        }
        i_domain_min = recvData[0];
        i_domain_max = recvData[1];
        j_domain_min = recvData[2];
        j_domain_max = recvData[3];
    }
}

void Communication::communicateNeighbourInfo(int my_rank, 
                                            int their_rank, 
                                            int &left_neighbour_rank,
                                            int &right_neighbour_rank, 
                                            int &bottom_neighbour_rank, 
                                            int &top_neighbour_rank)
{
    int mpi_status = MPI_SUCCESS;
    if (my_rank == 0){   
        int sendData[4] = {left_neighbour_rank, right_neighbour_rank, bottom_neighbour_rank,top_neighbour_rank};
        mpi_status = MPI_Send(&sendData, 4, MPI_INT, their_rank, 0, MPI_COMM_WORLD);
        if (MPI_SUCCESS != mpi_status){
            std::cerr << "MPI_Send fail.\n";
        }
    }
    else{
        int recvData[4] = {-1,-1,-1,-1};
        MPI_Status status;
        mpi_status = MPI_Recv(&recvData, 4, MPI_INT, their_rank, 0, MPI_COMM_WORLD, &status);
        if (MPI_SUCCESS != mpi_status){
            std::cerr << "MPI_Recv fail.\n";
        }
        left_neighbour_rank = recvData[0];
        right_neighbour_rank = recvData[1];
        bottom_neighbour_rank = recvData[2];
        top_neighbour_rank = recvData[3];
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

