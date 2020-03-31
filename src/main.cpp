#include <mpi.h>
#include <bits/stdc++.h>
#include <sys/resource.h>
#include <unistd.h>
#include "static/cloudmpi.h"

using namespace std;

#define MASTER 0

int main(int argc, char ** argv)
{
    int rank, ranks, size, times, window;
    int args[3];

    /* start up */
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);

    //TEST CODE
//    char processor_name[MPI_MAX_PROCESSOR_NAME];
//    int len;
//    MPI_Get_processor_name(processor_name, &len);
//    cout<<"Hello World ! from processor "<< processor_name<<" rank "<<rank<<" out of "<<ranks<<" processors\n";
//    cout<<"Sum :"<<add(5,6);
    //TEST CODE END

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
    if (rank == 0) {
        /* mark start of output */
        printf("START mpiGraph \n");
        printf("MsgSize\t%d\nTimes\t%d\nWindow\t%d\n",size,times,window);
        printf("Procs\t%d\n\n",ranks);
    }

    /* synchronize, then start the run */
    MPI_Barrier(MPI_COMM_WORLD);
    double ** arr = code(rank, ranks, size, times, window);
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0){
        printf("\n");
        printf("Combined\t\t\t");
        for(int k=0; k<ranks; k++) {
            //printf("%s:%d\t", &hostnames[k*sizeof(hostname)], k);
        }
        printf("\n");
        for(int j=0; j<ranks; j++) {
            //printf("%s:%d from\t\t", &hostnames[j*sizeof(hostname)], j);
            for(int k=0; k<ranks; k++) {
                printf("%0.3f\t\t", arr[j][k]);
            }
            printf("\n");
        }
    }
    /* print memory usage */

//  if(rank == 0) { printf("\n"); }
//  print_mpi_resources();
//

    /* mark end of output */
    if (rank == 0) { printf("END mpiGraph\n"); }

    // Test Rank allocation code
    MPI::Intracomm NodeComm;
    MPI::Intracomm MasterComm;

    int NodeRank, MasterRank, NodeSize, MasterSize;
    string NodeNameStr;
    bool b = CommByNode(NodeComm, MasterComm, NodeRank, MasterRank, NodeSize, MasterSize, NodeNameStr);

    //cout<<"\nBool value of rank "<<rank<<" : "<<b;
    cout<<'\n'<<NodeNameStr;
    printf("\nReal Rank : %d | Noderank : %d | NodeSize : %d",rank,NodeRank,NodeSize);
    printf("\nReal Rank : %d | Masterrank : %d | MasterSize : %d\n",rank,MasterRank,MasterSize);

    /* shut down */
    MPI_Finalize();
    return 0;
}