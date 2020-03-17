#include <mpi.h>
#include <bits/stdc++.h>
#include <sys/resource.h>
#include <unistd.h>
#include "cloudmpi.h"

using namespace std;

char  hostname[256];
char* hostnames;
char VERS[] = "1.5";


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

                int flag_sends = 0, flag_recvs = 0;

                while(!flag_sends || !flag_recvs) {
                    if (!flag_sends) {
                        MPI_Testall((k + 1) / 2, &request_array[(k + 1) / 2 - 1], &flag_sends,
                                    &status_array[(k + 1) / 2 - 1]);
                        if (flag_sends) {
                            TIME_END_SEND;
                        }
                    }
                    if (!flag_recvs) {
                        MPI_Testall((k + 1) / 2, &request_array[0], &flag_recvs, &status_array[(k + 1) / 2 - 1]);
                        if (flag_recvs) {
                            TIME_END_RECV;
                        }
                    }
                }
                sendtimes[send_rank*iter+i] = TIME_USECS_SEND / (double)w;
                recvtimes[recv_rank*iter+i] = TIME_USECS_RECV / (double)w;
                }

                dist++;

            }

        if(world_rank == 0)
            cout<<"Collecting Results\n";

        double *sendsums = (double*)malloc(sizeof(double)*world_size);
        double *recvsums = (double*)malloc(sizeof(double)*world_size);

        for(j=0;j<world_size;j++){
            sendsums[j] = 0.0;
            recvsums[j] = 0.0;

            if(j == world_rank)
                continue;

            for(i=0;i<iter;i++){
                double sendval = sendtimes[j*iter+i];
                double recvval = recvtimes[j*iter+i];

                sendsums[j] += sendval;
                recvsums[j] += recvval;
            }
        }

        double *allsums;

        if(world_rank == 0)
            allsums = (double*)malloc(sizeof(double)*world_size*world_size);

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(sendsums, world_rank, MPI_DOUBLE, allsums, world_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Send Stats
        if(world_rank == 0){
            double sendsum = 0.0;
            double sendmin = 10000000000000000.0;
            double MBpsec = ((double)(message_size))*1000000.0 / (1024.0 * 1024.0);

            for(j=0; j<world_size; j++){
                for(k=0; k<world_size; k++){
                    if(j==k)
                        continue;
                    double sendval = allsums[j*world_size+k];
                    sendsum += sendval;
                    sendmin = (sendval <sendmin) ? sendval : sendmin;
                }
            }

            //Printing Sent Stats

            sendmin /= (double) iter;
            sendsum /= (double) (world_size)*(world_size-1)*iter;

            cout<<"\nSend max : "<<MBpsec/sendmin;
            cout<<"\nSend avg : "<<MBpsec/sendsum;

            //Printing Bandwidth Table

            cout<<"\nSend\t";

            for(k=0;k<world_size;k++){
                printf("%s:%d\t", &hostnames[k*sizeof(hostname)], k);
            }
            cout<<"\n";

            for(j=0;j<world_size;j++){
                printf("%s:%d to \t", &hostnames[j*sizeof(hostname)], j);

                for(k=0;k<world_size;k++){
                    double val = allsums[j*world_size+k];

                    if(val != 0.0){
                        val = MBpsec * (double) iter /val ;
                    }
                    printf("%0.3f\t", val);
                }
                cout<<"\n";
            }
        }


        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(recvsums, world_size, MPI_DOUBLE, allsums, world_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        //Recv Stats
        if(world_rank == 0){
            double recvsum =0.0;
            double recvmin = 10000000000000000.0;
            double MBpsec = ((double)(message_size))*1000000.0 / (1024.0*1024.0);

            for(j=0;j<world_size;j++){
                for(k=0;k<world_size;k++){
                    if(j == k)
                        continue;
                    double recvval = allsums[j*world_size+k];
                    recvsum += recvval;
                    recvmin = (recvval < recvmin) ? recvval : recvmin;
                }
            }

            //Print recv stats
            recvmin /=(double) iter;
            recvsum /= (double) (world_size)*(world_size-1)*iter;
            cout<<"\nRecv max : "<<MBpsec/recvmin;
            cout<<"\nRecv avg : "<<MBpsec/recvsum;

            //Print recv bandwidth
            cout<<"\nRecv\t";

            for(k=0;k<world_size;k++){
                printf("%s : %d\t", &hostnames[k* sizeof(hostname)], k);
            }
            cout<<"\n";
            for(j=0;j<world_size;j++){
                printf("%s:%d", &hostnames[k*sizeof(hostname)], k);
                for(k=0;k<world_size;k++){
                    double val = allsums[j*world_size+k];
                    if(val != 0.0){
                        val = MBpsec* (double) iter/val;
                    }
                    printf("%0.3f\t", val);
                }
                cout<<"\n";
            }
        }

        //Free of memory
        if(world_rank == 0){
            free(allsums);
        }
        free(sendsums);
        free(recvsums);
        free(send);
        free(recv);
        free(status_array);
        free(request_array);
        free(sendtimes);
        free(recvtimes);
        free(message_tags);

        return;
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

    gethostname(hostname, sizeof(hostname));
    hostnames = (char*)malloc(sizeof(hostname)*world_size);
    MPI_Gather(hostname, sizeof(hostname), MPI_CHAR, hostnames, sizeof(hostname), MPI_CHAR,0,MPI_COMM_WORLD);

    int message_size, iter, window, args[3];
    message_size = 4096*4;
    iter = 10;
    window = 50;

    if(argc == 4){
        message_size = atoi(argv[1]);
        iter = atoi(argv[2]);
        window = atoi(argv[3]);
    }
    args[0] = message_size;
    args[1] = iter;
    args[2] = window;

    if(world_rank == 0){
        cout<<"Start mpiGraph "<<VERS;
        cout<<"\nMsgSize : "<<message_size<<"\nTimes : "<<iter<<"\nWindow : "<<window;
        cout<<"\nTotal Processes : "<<world_size;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    pingpong(world_rank,world_size,message_size,iter,window);

    if(world_rank == 0){
        cout<<"Process Complete";
    }

    MPI_Finalize();
    return 0;

        

    MPI_Finalize();
}