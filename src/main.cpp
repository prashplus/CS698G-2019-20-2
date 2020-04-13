#include <mpi.h>
#include <bits/stdc++.h>
#include <sys/resource.h>
#include <unistd.h>
#include "static/cloudmpi.h"

using namespace std;

#define MASTER 0

int main(int argc, char ** argv)
{
    int world_rank, world_size, size, times, window;
    int args[3];

    /* start up */
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /* set job parameters, read values from command line if they're there */
    size = 1000000;
    times = 100;
    window = 50;
    if (argc == 4) {
        size   = atoi(argv[1]);
        times  = atoi(argv[2]);
        window = atoi(argv[3]);
    }
    args[0] = size;
    args[1] = times;
    args[2] = window;

    /* print the header */
    if (world_rank == 0) {
        /* mark start of output */
        printf("START mpiGraph \n");
        printf("MsgSize\t%d\nTimes\t%d\nWindow\t%d\n",size,times,window);
        printf("Procs\t%d\n\n",world_size);
    }

    /* synchronize, then start the run */
    MPI_Barrier(MPI_COMM_WORLD);
    double ** dist = code(world_rank, world_size, size, times, window);
    MPI_Barrier(MPI_COMM_WORLD);

    if(world_rank == 0){
        printf("\n");
        printf("Combined\t\t\t");
        for(int k=0; k<world_size; k++) {
            //printf("%s:%d\t", &hostnames[k*sizeof(hostname)], k);
        }
        printf("\n");
        for(int j=0; j<world_size; j++) {
            //printf("%s:%d from\t\t", &hostnames[j*sizeof(hostname)], j);
            for(int k=0; k<world_size; k++) {
                printf("%0.3f\t\t", dist[j][k]);
            }
            printf("\n");
        }
    }

    /* mark end of output */
    if (world_rank == 0) { printf("END mpiGraph\n"); }

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
//    if(l2_NodeRank == 0)
//        printf("\nMYBCAST: Real Rank : %d | NodeRank : %d | MasterRank : %d | Root Rank : %d | Time : %lf", rank, NodeRank, MasterRank, root_rank, time1);

    MPI_Barrier(MPI_COMM_WORLD);
    time2 -= MPI_Wtime();
    MPI_Bcast(data, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    time2 += MPI_Wtime();
//    if(world_rank == 0)
//        printf("\nBCAST: Real Rank : %d | NodeRank : %d | MasterRank : %d | Root Rank : %d | Time : %lf\n", rank, NodeRank, MasterRank, root_rank, time2);



    l1_create_comm(l1_NodeComm,l1_MasterComm, l1_root_comm, dist);
    free(data);
    /* shut down */

    //TEST CODE
//    char processor_name[MPI_MAX_PROCESSOR_NAME];
//    int len;
//    MPI_Get_processor_name(processor_name, &len);
//    cout<<"Hello World ! from processor "<< processor_name<<" rank "<<rank<<" out of "<<ranks<<" processors\n";
//    cout<<"Sum :"<<add(5,6);
    //TEST CODE END

    MPI_Finalize();
    return 0;
}