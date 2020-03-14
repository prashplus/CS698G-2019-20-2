#include <mpi.h>
#include <bits/stdc++.h>
#include "cloudmpi.h"
using namespace std;


//TIming calculations
#define TIME_START (time_start = MPI_Wtime())
#define TIME_END_SEND (time_end_send = MPI_Wtime())
#define TIME_END_RECV (time_end_recv = MPI_Wtime())
#define TIME_USECS_SEND ((time_end_send - time_start) * 1000000.0)
#define TIME_USECS_RECV ((time_end_recv - time_start) * 1000000.0)

double time_start, time_end_send, time_end_recv;

void pingpong(int world_rank, int world_size, int message_size, int iter, int window)
    {                                                                                                                                                                        
        int i,j,k,w;

        char* send = (char*)malloc(window*message_size);
        char* recv = (char*)malloc(window*message_size);

        MPI_Status*  status_array  = (MPI_Status*)  malloc(sizeof(MPI_Status) *window*2);
        MPI_Request* request_array = (MPI_Request*) malloc(sizeof(MPI_Request)*window*2);
        double* sendtimes = (double*) malloc(sizeof(double)*iter*world_size);
        double* recvtimes = (double*) malloc(sizeof(double)*iter*world_size);

        int* message_tags = (int*)malloc(window* sizeof(int));
        for (i=0;j<window;i++){
            message_tags[i]=i;
        }

        int dist = 1;
        while(dist < world_size) {
            float progress = (float) dist/ (float) world_size * 100.0;
            if(world_rank == 0) {
                printf("%d of %d (%0.1f%%)\n", dist, world_size, progress);
                fflush(stdout);
            }

            int send_rank = (world_rank + dist + world_size) % world_size;
            int recv_rank = (world_rank - dist + world_size) % world_size;

            for(i=0;i<iter;i++){
                //Sync steps
                MPI_Barrier(MPI_COMM_WORLD);

                TIME_START;
                k = -1;

                for(w=0; w<window; w++){
                    k++;
                    MPI_Irecv(&recv[w*message_size], message_size, MPI_BYTE, recv_rank, MPI_ANY_TAG, MPI_COMM_WORLD, &request_array[k]);
                }

                for(w=0; w<window; w++){
                    k++;
                    MPI_Isend(&send[w*message_size], message_size, MPI_BYTE, recv_rank, MPI_ANY_TAG, MPI_COMM_WORLD, &request_array[k]);
                }

                int flag_sends = 0, flag_recvs =0;

//                while(!flag_sends || !flag_recvs){
//
//                }

            }

        }


    }

int main(int argc, char ** argv)
{
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //TEST CODE
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name(processor_name, &len);
    cout<<"Hello World ! from processor "<< processor_name<<" rank "<<world_rank<<" out of "<<world_size<<" processors\n";
    cout<<"Sum :"<<add(5,6);
    //TEST CODE END

        

    MPI_Finalize();
}