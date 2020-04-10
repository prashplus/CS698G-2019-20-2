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
//    if (rank == 0) {
//        /* mark start of output */
//        printf("START mpiGraph \n");
//        printf("MsgSize\t%d\nTimes\t%d\nWindow\t%d\n",size,times,window);
//        printf("Procs\t%d\n\n",ranks);
//    }
//
//    /* synchronize, then start the run */
//    MPI_Barrier(MPI_COMM_WORLD);
//    double ** arr = code(rank, ranks, size, times, window);
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    if(rank == 0){
//        printf("\n");
//        printf("Combined\t\t\t");
//        for(int k=0; k<ranks; k++) {
//            //printf("%s:%d\t", &hostnames[k*sizeof(hostname)], k);
//        }
//        printf("\n");
//        for(int j=0; j<ranks; j++) {
//            //printf("%s:%d from\t\t", &hostnames[j*sizeof(hostname)], j);
//            for(int k=0; k<ranks; k++) {
//                printf("%0.3f\t\t", arr[j][k]);
//            }
//            printf("\n");
//        }
//    }
//
//    /* mark end of output */
//    if (rank == 0) { printf("END mpiGraph\n"); }
//
    // Test Rank allocation code
    MPI::Intracomm NodeComm;
    MPI::Intracomm MasterComm;

    int NodeRank, MasterRank, NodeSize, MasterSize;
    string NodeNameStr;
    bool b = CommByNode(NodeComm, MasterComm, NodeRank, MasterRank, NodeSize, MasterSize, NodeNameStr);

//    //cout<<"\nBool value of rank "<<rank<<" : "<<b;
//    cout<<'\n'<<NodeNameStr;
//    printf("\nReal Rank : %d | NodeRank : %d | MasterRank : %d",rank,NodeRank,MasterRank);
//    printf("\nReal Size : %d | NodeSize : %d | MasterSize : %d\n",ranks,NodeSize,MasterSize);
//    //Matrix end

    // Get the group of processes in MPI_COMM_WORLD
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);

    int r[MasterSize],i,j,k;
    int *rankMatrix,*rankData;

//    if(rank == 0){
//        rankMatrix = (int **)malloc(ranks * sizeof(int *));
//        for (i=0; i<ranks; i++)
//            rankMatrix[i] = (int *)malloc(3 * sizeof(int));
//    }
    rankMatrix = (int*)malloc(sizeof(int)*3*ranks);
    rankData = (int*)malloc(sizeof(int)*3);
    rankData[0]=rank;
    rankData[1]=NodeRank;
    rankData[2]=MasterRank;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(rankData, 3, MPI_INT, rankMatrix, 3, MPI_INT, 0, MPI_COMM_WORLD);

    if(rank == 0){
        printf("\n");
        printf("\t\tRank Data\n");
        printf("\n");
        printf("World Rank\tNode Rank\tMaster Rank\n");
        for(j=0; j<ranks; j++) {
            for(k=0; k<3; k++) {
                int val = rankMatrix[j*3+k];
                printf("%d\t\t", val);
            }
            printf("\n");
        }
    }

    /* shut down */
    MPI_Finalize();
    return 0;
}