#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char ** argv)
{
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name(processor_name, &len);

    cout<<"Hello World ! from processor "<< processor_name<<" rank "<<world_rank<<" out of "<<world_size<<" processors\n";

    MPI_Finalize();
}