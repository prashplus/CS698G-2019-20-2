#include <mpi.h>
#include <bits/stdc++.h>
#include <sys/resource.h>
#include <unistd.h>
#include "static/cloudmpi.h"

using namespace std;

#define MASTER 0

long size, times, window;
int world_rank, world_size;
double THRESHOLD;

void test1(){
    //TEST CODE
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name(processor_name, &len);
    cout<<"Hello World ! from processor "<< processor_name<<" rank "<<world_rank<<" out of "<<world_size<<" processors\n";
    cout<<"Sum :"<<add(5,6);
    //TEST CODE END
}

double ** test2(){
    /* print the header */
    if (world_rank == 0) {
        /* mark start of output */
        printf("START mpiGraph \n");
        printf("MsgSize\t%ld\nTimes\t%ld\nWindow\t%ld\n",size,times,window);
        printf("Procs\t%d\n\n",world_size);
    }

    /* synchronize, then start the run */
    MPI_Barrier(MPI_COMM_WORLD);
    double ** dist = code(world_rank, world_size, size, times, window);
    MPI_Barrier(MPI_COMM_WORLD);

    /* mark end of output */
    if (world_rank == 0) { printf("END mpiGraph\n"); }
    return dist;
}

int main(int argc, char ** argv)
{

    /* start up */
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    char *end;
    /* set job parameters, read values from command line if they're there */
    size = 10000;
    times = 10;
    window = 10;
    THRESHOLD = 500.0;

    if (argc == 5) {
        size   = strtol(argv[1],&end,10);
        times  = strtol(argv[2],&end,10);
        window = strtol(argv[3],&end,10);
        THRESHOLD = strtol(argv[4],&end,10);
    }

    //Test 1 : Processor Name
    //test1();

    //Test 2 : Check Latency Matrix
    double **dist = test2();

    // Test Rank allocation code
    MPI::Intracomm l1_NodeComm;
    MPI::Intracomm l1_MasterComm;
    MPI_Comm l1_root_comm;
    MPI::Intracomm l2_NodeComm;
    MPI::Intracomm l2_MasterComm;
    MPI_Comm l2_root_comm;


    l2_create_comm(l2_NodeComm,l2_MasterComm, l2_root_comm);
    double time1,time2,*data;
    data = (double *)malloc(sizeof(double)*size);
    time1 -= MPI_Wtime();
    if(MPI_COMM_NULL != l2_root_comm){
        MPI_Barrier(l2_root_comm);
        MPI_Bcast(data, size, MPI_DOUBLE, 0, l2_root_comm);
        MPI_Barrier(l2_root_comm);
    }

    MPI_Barrier(l2_NodeComm);
    MPI_Bcast(data, size, MPI_DOUBLE, 0, l2_NodeComm);
    MPI_Barrier(l2_NodeComm);
    time1 += MPI_Wtime();
    if(world_rank == 0)
        printf("\nCustomBCAST: Real Rank : %d | Time : %lf", world_rank, time1);

    MPI_Barrier(MPI_COMM_WORLD);
    time2 -= MPI_Wtime();
    MPI_Bcast(data, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    time2 += MPI_Wtime();
    if(world_rank == 0)
        printf("\nMPI_BCAST: Real Rank : %d | Time : %lf", world_rank, time2);
    l1_create_comm(l1_NodeComm,l1_MasterComm, l1_root_comm, dist, THRESHOLD);

    free(data);
    /* shut down */


    MPI_Finalize();
    return 0;
}