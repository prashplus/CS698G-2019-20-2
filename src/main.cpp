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
    if(world_rank == 0){
        printf("\n");
        printf("Combined\t\t\t");
//        for(int k=0; k<world_size; k++) {
//            printf("%s:%d\t", &hostnames[k*sizeof(hostname)], k);
//        }
        printf("\n");
        for(int j=0; j<world_size; j++) {
//            printf("%s:%d from\t\t", &hostnames[j*sizeof(hostname)], j);
            for(int k=0; k<world_size; k++) {
                printf("%0.3f\t\t", dist[j][k]);
            }
            printf("\n");
        }
    }
}

void test3(int size){
    double time1=0,time2=0,ftime1,ftime2,*data;
    data = (double *)malloc(sizeof(double)*size);

    // Init part of Data at root
    if(world_rank == 0){
        for(int i = 0;i < 4; i++)
            data[i]=i;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    time1 -= MPI_Wtime();
    MPI_CustomBcast(data, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
    time1 += MPI_Wtime();
    MPI_Reduce(&time1, &ftime1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

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
    MPI_Reduce(&time2, &ftime2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(world_rank == 0)
    {
        printf("\nMPI_CustomBcast Time : %lf", ftime1);
        printf("\nMPI_BCAST Time : %lf", ftime2);

        if (!ifstream("data.txt"))
        {
            ofstream temp;
            temp.open("data.txt",ios::out);
            temp.close();
        }
        ofstream f;
        f.open ("data.txt",ios::app);
        f<<ftime1<<'\n';
        f<<ftime2<<'\n';
        f.close();
    }

    free(data);
}

void test4(int size){
    double time1=0,time2=0,ftime1,ftime2,*send_data, *recv_data;
    int send_count = size,recv_count = size;
    if(world_rank == 0)
        send_data = (double *)malloc(sizeof(double)*send_count);

    recv_data = (double *)malloc(sizeof(double)*recv_count);
    // Init part of Data at root
    if(world_rank == 0){
        for(int i = 0;i < 4; i++)
            send_data[i]=i;
    }

    //Custom Scatter
    MPI_Barrier(MPI_COMM_WORLD);
    time1 -= MPI_Wtime();
    MPI_CustomScatter(send_data, send_count, MPI_DOUBLE, recv_data, recv_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    time1 += MPI_Wtime();
    MPI_Reduce(&time1, &ftime1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    cout<<"\nRank : "<<world_rank<<"\n ";
    cout<<recv_data[0]<<"\t";

    //BUILT IN Scatter
    MPI_Barrier(MPI_COMM_WORLD);
    time2 -= MPI_Wtime();
    MPI_Scatter(send_data, send_count, MPI_DOUBLE, recv_data, recv_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    time2 += MPI_Wtime();
    MPI_Reduce(&time2, &ftime2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(world_rank == 0)
    {
        printf("\nMPI_CustomScatter Time : %lf", ftime1);
        printf("\nMPI_Scatter Time : %lf", ftime2);

        if (!ifstream("scatter.txt"))
        {
            ofstream temp;
            temp.open("scatter.txt",ios::out);
            temp.close();
        }
        ofstream f;
        f.open ("scatter.txt",ios::app);
        f<<ftime1<<'\n';
        f<<ftime2<<'\n';
        f.close();
    }

    free(send_data);
    free(recv_data);
}

void test5(int size){
    double time1=0,time2=0,ftime1,ftime2;
    int *send_data, *recv_data;

    send_data = (int *)malloc(sizeof(int)*size);
    recv_data = (int *)malloc(sizeof(int)*size);

    for(int i = 0;i < size;i++){
        send_data[i]=1;
    }

    //Custom Allreduce
    MPI_Barrier(MPI_COMM_WORLD);
    time1 -= MPI_Wtime();
    MPI_CustomAllreduce(send_data, recv_data, size, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    time1 += MPI_Wtime();
    MPI_Reduce(&time1, &ftime1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    cout<<"\nRank : "<<world_rank<<"\n ";
    for(int i=0;i<5;i++)
        cout<<recv_data[i]<<"\t";

    MPI_Barrier(MPI_COMM_WORLD);
    time2 -= MPI_Wtime();
    MPI_Allreduce(send_data, recv_data, size, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    time2 += MPI_Wtime();
    MPI_Reduce(&time2, &ftime2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


    if(world_rank == 0)
    {
        printf("\nMPI_CustomAllreduce Time : %lf", ftime1);
        printf("\nMPI_Allreduce Time : %lf", ftime2);

        if (!ifstream("allreduce.txt"))
        {
            ofstream temp;
            temp.open("allreduced.txt",ios::out);
            temp.close();
        }
        ofstream f;
        f.open ("allreduce.txt",ios::app);
        f<<ftime1<<'\n';
        f<<ftime2<<'\n';
        f.close();
    }

//    free(send_data);
//    free(recv_data);
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
        printf("Initializing (will take couple of minutes)\n");
        printf("Parameters :\nProcesses : \t%d\nData Size : \t%ld\n",world_size,size);
    }


    /* **************************** YOU MUST CALL init() to use any of Custom Function *****************************
     *
     *  */
    //Initialize all the Comms
    init();


    /* ********************************** All the test cases of Custom Functions ***********************************
     *
     *  */
//    //Test 1 : Processor Name
//    test1();

//    //Test 2 : Check Latency Matrix
//    test2();

    //Test 3 : MPI_CustomBcast
    test3(128);
    test3(1280);
    test3(12800);
    test3(131072);
    test3(1310720);
//    test3(10000000);
//    test3(100000000);
//    test3(1000000000);

//    //Test 4: MPI_CustomScatter
//    test4(size);

//    //Test 5: MPI_CustomALlreduce
//    test5(size);

    MPI_Finalize();
    return 0;
}