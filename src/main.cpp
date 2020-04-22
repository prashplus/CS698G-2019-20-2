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

void test2(){
    double **dist = getDist();
}

void test3(int size){
    double time1,time2,*data;
    data = (double *)malloc(sizeof(double)*size);

    // Init part of Data at root
    if(world_rank == 0){
        for(int i = 0;i < 4; i++)
            data[i]=i;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    time1 -= MPI_Wtime();
    MPI_CustomBcast(data, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    time1 += MPI_Wtime();
    if(world_rank == 0)
        printf("\nMPI_CustomBcast: Real Rank : %d | Time : %lf", world_rank, time1);

    // Start Data Received
//    printf("\n Data Received at rank : %d | ",world_rank);
//    for(int i =0;i<4;i++){
//        printf(" %lf ",data[i]);
//    }
//    printf("\n");
    // End Data received


    // BuiltIn Bcast
    MPI_Barrier(MPI_COMM_WORLD);
    time2 -= MPI_Wtime();
    MPI_Bcast(data, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    time2 += MPI_Wtime();
    if(world_rank == 0)
        printf("\nMPI_BCAST: Real Rank : %d | Time : %lf", world_rank, time2);

    free(data);
}

int main(int argc, char ** argv)
{
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    char *end;

    // Command line parameters
    size = 1000;
    times = 10;
    window = 10;
    THRESHOLD = 500.0;

    if (argc == 2) {
        size   = strtol(argv[1],&end,10);
//        times  = strtol(argv[2],&end,10);
//        window = strtol(argv[3],&end,10);
//        THRESHOLD = strtol(argv[4],&end,10);
    }
    if (world_rank == 0) {
        printf("START mpiGraph \n");
        printf("MsgSize\t%ld\nTimes\t%ld\nWindow\t%ld\n",size,times,window);
        printf("Procs\t%d\n\n",world_size);
    }
    //Initialize all the Comms
    init();
//    //Test 1 : Processor Name
//    test1();
//
//    //Test 2 : Check Latency Matrix
//    test2();
//
//    //Test 3 : MPI_CustomBcast
    test3(size);
    /* shut down */

    MPI_Finalize();
    return 0;
}