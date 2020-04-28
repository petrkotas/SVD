// Include the MPI version 2 C++ bindings:
#include <mpi.h>
#include <iostream>
#include <string.h>

using namespace std;

int rank, pr; 

int main(int argc, char* argv[])
{
    /*Inicializacni rutina pro zavedeni MPI*/
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &pr);

    cout << "hello_parallel.cc: Number of tasks="<<pr<<" My rank=" << rank <<endl;

    // Tell the MPI library to release all resources it is using:
    MPI_Finalize();
    return 0;
}
