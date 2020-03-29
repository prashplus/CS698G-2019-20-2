#include <mpi.h>
#include <bits/stdc++.h>
#include <sys/resource.h>
#include <unistd.h>
#include "static/cloudmpi.h"

using namespace std;

char  hostname[256];
char* hostnames;
char VERS[] = "1.5";
#define MASTER 0


/* =============================================================
 * TIMER MACROS
 * These macros start/stop the timer and measure the difference
 * =============================================================
 */

#ifdef USE_GETTIMEOFDAY
/* use gettimeofday() for timers */

#include <sys/time.h>
#define __TIME_START__    (gettimeofday(&g_timeval__start, &g_timezone))
#define __TIME_END_SEND__ (gettimeofday(&g_timeval__end_send, &g_timezone))
#define __TIME_END_RECV__ (gettimeofday(&g_timeval__end_recv, &g_timezone))
#define __TIME_USECS_SEND__ (d_Time_Diff_Micros(g_timeval__start, g_timeval__end))
#define __TIME_USECS_RECV__ (d_Time_Diff_Micros(g_timeval__start, g_timeval__end_recv))
#define d_Time_Diff_Micros(timeval__start, timeval__end) \
  ( \
    (double) (  (timeval__end.tv_sec  - timeval__start.tv_sec ) * 1000000 \
              + (timeval__end.tv_usec - timeval__start.tv_usec)  ) \
  )
#define d_Time_Micros(timeval) \
  ( \
    (double) (  timeval.tv_sec * 1000000 \
              + timeval.tv_usec  ) \
  )
struct timeval  g_timeval__start, g_timeval__end_send, g_timeval__end_recv;
struct timezone g_timezone;

#else
/* use MPI_Wtime() for timers instead of gettimeofday() (recommended)
 * on some systems gettimeofday may reset backwards via some global clock,
 * which leads to incorrect timing data including negative time periods
 */

#define __TIME_START__    (g_timeval__start    = MPI_Wtime())
#define __TIME_END_SEND__ (g_timeval__end_send = MPI_Wtime())
#define __TIME_END_RECV__ (g_timeval__end_recv = MPI_Wtime())
#define __TIME_USECS_SEND__ ((g_timeval__end_send - g_timeval__start) * 1000000.0)
#define __TIME_USECS_RECV__ ((g_timeval__end_recv - g_timeval__start) * 1000000.0)
double g_timeval__start, g_timeval__end_send, g_timeval__end_recv;

#endif /* of USE_GETTIMEOFDAY */


/* =============================================================
 * MAIN TIMING LOGIC
 * Uses a ring-based (aka. shift-based) algorithm.
 * 1) First, logically arrange the MPI tasks
 *    from left to right from rank 0 to rank N-1 in a circular array.
 * 2) Then, during each step, each task sends messages to the task D uints to the right
 *    and receives from task D units to the left.
 *    In each step, each tasks measures its send and receive bandwidths.
 * 3) There are N-1 such steps so that each task has sent to and received from every task.
 * =============================================================
 */
void code(int mypid, int nnodes, int size, int times, int window)
{
    /* arguments are:
     *   mypid  = rank of this process
     *   nnodes = number of ranks
     *   size   = message size in bytes
     *   times  = number of times to measure bandwidth between task pairs
     *   window = number of outstanding sends and recvs to a single rank
     */
    int i, j, k, w;

    /* allocate memory for all of the messages */
    char* send_message = (char*) malloc(window*size);
    char* recv_message = (char*) malloc(window*size);
    MPI_Status*  status_array  = (MPI_Status*)  malloc(sizeof(MPI_Status) *window*2);
    MPI_Request* request_array = (MPI_Request*) malloc(sizeof(MPI_Request)*window*2);
    double* sendtimes = (double*) malloc(sizeof(double)*times*nnodes);
    double* recvtimes = (double*) malloc(sizeof(double)*times*nnodes);

    int* message_tags = (int*) malloc(window*sizeof(int));
    for (i=0;i<window;i++) { message_tags[i] = i; }

    /* start iterating over distance */
    int distance = 1;
    while (distance < nnodes) {
        /* this test can run for a long time, so print progress to screen as we go */
        float progress = (float) distance / (float) nnodes * 100.0;
        if (mypid == 0) {
            printf("%d of %d (%0.1f%%)\n", distance, nnodes, progress);
            fflush(stdout);
        }

        /* find tasks distance units to the right (send) and left (recv) */
        int sendpid = (mypid + distance + nnodes) % nnodes;
        int recvpid = (mypid - distance + nnodes) % nnodes;

        /* run through 'times' iterations on a given ring */
        for (i=0; i<times; i++) {
            /* couple of synch's to make sure everyone is ready to go */
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            __TIME_START__;
            k=-1;
            /* setup a window of irecvs from my partner who is distance steps to my left */
            for (w=0; w<window; w++) {
                k=k+1;
                MPI_Irecv(&recv_message[w*size], size, MPI_BYTE,
                          recvpid, MPI_ANY_TAG, MPI_COMM_WORLD, &request_array[k]);
            }
            /* fire off a window of isends to my send partner distance steps to my right */
            for (w=0; w<window; w++) {
                k=k+1;
                MPI_Isend(&send_message[w*size], size, MPI_BYTE,
                          sendpid, message_tags[w], MPI_COMM_WORLD, &request_array[k]);
            }
            /* time sends and receives separately */
            int flag_sends = 0;
            int flag_recvs = 0;
            while(!flag_sends || !flag_recvs) {
                /* check whether the sends are done */
                if (!flag_sends) {
                    MPI_Testall((k+1)/2, &request_array[(k+1)/2-1], &flag_sends, &status_array[(k+1)/2-1]);
                    if (flag_sends) { __TIME_END_SEND__; }
                }

                /* check whether the recvs are done */
                if (!flag_recvs) {
                    MPI_Testall((k+1)/2, &request_array[0], &flag_recvs, &status_array[0]);
                    if (flag_recvs) { __TIME_END_RECV__; }
                }
            }
            sendtimes[sendpid*times+i] = __TIME_USECS_SEND__ / (double) w;
            recvtimes[recvpid*times+i] = __TIME_USECS_RECV__ / (double) w;
        } /* end times loop */
        /* bump up the distance for the next ring step */
        distance++;
    } /* end distance loop */

    /* for each node, compute sum of my bandwidths with that node */
    if(mypid == 0) printf("Gathering results\n");
    double* sendsums = (double*) malloc(sizeof(double)*nnodes);
    double* recvsums = (double*) malloc(sizeof(double)*nnodes);
    for(j=0; j<nnodes; j++) {
        sendsums[j] = 0.0;
        recvsums[j] = 0.0;
        if (j == mypid) continue;
        for(i=0; i<times; i++) {
            double sendval = sendtimes[j*times+i];
            double recvval = recvtimes[j*times+i];
            sendsums[j] += sendval;
            recvsums[j] += recvval;
        }
    }

    /* gather send bw sums to rank 0 */
    double* allsums;
    if (mypid == 0) {
        allsums = (double*) malloc(sizeof(double)*nnodes*nnodes);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(sendsums, nnodes, MPI_DOUBLE, allsums, nnodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int r = nnodes;
    double **arr = (double **)malloc(r * sizeof(double *));
    for (i=0; i<r; i++)
        arr[i] = (double *)malloc(r * sizeof(double));
    /* rank 0 computes send stats and prints result */
    if (mypid == 0) {
        /* compute stats over all nodes */
        double sendsum = 0.0;
        double sendmin = 10000000000000000.0;
        double MBsec   = ((double)(size)) * 1000000.0 / (1024.0*1024.0);
        for(j=0; j<nnodes; j++) {
            for(k=0; k<nnodes; k++) {
                if (j == k) continue;
                double sendval = allsums[j*nnodes+k];
                sendsum += sendval;
                sendmin = (sendval < sendmin) ? sendval : sendmin;
            }
        }

        /* print send stats */
        sendmin /= (double) times;
        sendsum /= (double) (nnodes)*(nnodes-1)*times;
        printf("\nSend max\t%f\n", MBsec/sendmin);
        printf("Send avg\t%f\n", MBsec/sendsum);

        /* print send bandwidth table */
        printf("\n");
        printf("Send\t\t\t");
        for(k=0; k<nnodes; k++) {
            printf("%s:%d\t", &hostnames[k*sizeof(hostname)], k);
        }
        printf("\n");
        for(j=0; j<nnodes; j++) {
            printf("%s:%d to\t\t", &hostnames[j*sizeof(hostname)], j);
            for(k=0; k<nnodes; k++) {
                double val = allsums[j*nnodes+k];
                if (val != 0.0) { val = MBsec * (double) times / val; }
                arr[j][k]=val;
                printf("%0.3f\t\t", val);
            }
            printf("\n");
        }
    }

    /* gather recv bw sums to rank 0 */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(recvsums, nnodes, MPI_DOUBLE, allsums, nnodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* rank 0 computes recv stats and prints result */
    if (mypid == 0) {
        /* compute stats over all nodes */
        double recvsum = 0.0;
        double recvmin = 10000000000000000.0;
        double MBsec   = ((double)(size)) * 1000000.0 / (1024.0*1024.0);
        for(j=0; j<nnodes; j++) {
            for(k=0; k<nnodes; k++) {
                if (j == k) continue;
                double recvval = allsums[j*nnodes+k];
                recvsum += recvval;
                recvmin = (recvval < recvmin) ? recvval : recvmin;
            }
        }

        /* print receive stats */
        recvmin /= (double) times;
        recvsum /= (double) (nnodes)*(nnodes-1)*times;
        printf("\nRecv max\t%f\n", MBsec/recvmin);
        printf("Recv avg\t%f\n", MBsec/recvsum);

        /* print receive bandwidth table */
        printf("\n");
        printf("Recv\t\t\t");
        for(k=0; k<nnodes; k++) {
            printf("%s:%d\t", &hostnames[k*sizeof(hostname)], k);
        }
        printf("\n");
        for(j=0; j<nnodes; j++) {
            printf("%s:%d from\t\t", &hostnames[j*sizeof(hostname)], j);
            for(k=0; k<nnodes; k++) {
                double val = allsums[j*nnodes+k];
                if (val != 0.0) { val = MBsec * (double) times / val; }
                arr[j][k]+=val;
                printf("%0.3f\t\t", val);
            }
            printf("\n");
        }

        /* Combined Matrix*/
        printf("\n");
        printf("Recv\t\t\t");
        for(k=0; k<nnodes; k++) {
            printf("%s:%d\t", &hostnames[k*sizeof(hostname)], k);
        }
        printf("\n");
        for(j=0; j<nnodes; j++) {
            printf("%s:%d from\t\t", &hostnames[j*sizeof(hostname)], j);
            for(k=0; k<nnodes; k++) {
                printf("%0.3f\t\t", arr[j][k]);
            }
            printf("\n");
        }
    }


    /* free off memory */
    if (mypid == 0) {
        free(allsums);
    }
    free(sendsums);
    free(recvsums);
    free(send_message);
    free(recv_message);
    free(status_array);
    free(request_array);
    free(sendtimes);
    free(recvtimes);
    free(message_tags);

    return;
}

bool CommByNode(MPI::Intracomm &NodeComm, MPI::Intracomm &MasterComm,
                   int &NodeRank, int &MasterRank, int &NodeSize, int &MasterSize,
                   string &NodeNameStr)
{
    bool IsOk = true;

    int Rank = MPI::COMM_WORLD.Get_rank();
    int Size = MPI::COMM_WORLD.Get_size();

    /*
     * ======================================================================
     * What follows is my best attempt at creating a communicator
     * for each node in a job such that only the cores on that
     * node are in the node's communicator, and each core groups
     * itself and the node communicator is made using the Split() function.
     * The end of this (lengthly) process is indicated by another comment.
     * ======================================================================
     */
    char *NodeName, *NodeNameList;
    NodeName = new char [1000];
    int NodeNameLen,
            *NodeNameCountVect,
            *NodeNameOffsetVect,
            NodeNameTotalLen = 0;
    //  Get the name and name character count of each core's node
    MPI::Get_processor_name(NodeName, NodeNameLen);

    //  Prepare a vector for character counts of node names
    if (Rank == MASTER)
        NodeNameCountVect = new int [Size];

    //  Gather node name lengths to master to prepare c-array
    MPI::COMM_WORLD.Gather(&NodeNameLen, 1, MPI::INT, NodeNameCountVect, 1, MPI::INT, MASTER);

    if (Rank == MASTER){
        //  Need character count information for navigating node name c-array
        NodeNameOffsetVect = new int [Size];
        NodeNameOffsetVect[0] = 0;
        NodeNameTotalLen = NodeNameCountVect[0];

        //  build offset vector and total char count for all node names
        for (int i = 1 ; i < Size ; ++i){
            NodeNameOffsetVect[i] = NodeNameCountVect[i-1] + NodeNameOffsetVect[i-1];
            NodeNameTotalLen += NodeNameCountVect[i];
        }
        //  char-array for all node names
        NodeNameList = new char [NodeNameTotalLen];
    }

    //  Gatherv node names to char-array in master
    MPI::COMM_WORLD.Gatherv(NodeName, NodeNameLen, MPI::CHAR, NodeNameList, NodeNameCountVect, NodeNameOffsetVect, MPI::CHAR, MASTER);

    string *FullStrList, *NodeStrList;
    //  Each core keeps its node's name in a str for later comparison
    stringstream ss;
    ss << NodeName;
    ss >> NodeNameStr;

    delete NodeName;    //  node name in str, so delete c-array

    int *NodeListLenVect, NumUniqueNodes = 0, NodeListCharLen = 0;
    string NodeListStr;

    if (Rank == MASTER){
        /*
         * Need to prepare a list of all unique node names, so first
         * need all node names (incl duplicates) as strings, then
         * can make a list of all unique node names.
         */
        FullStrList = new string [Size];    //  full list of node names, each will be checked
        NodeStrList = new string [Size];    //  list of unique node names, used for checking above list
        //  i loops over node names, j loops over characters for each node name.
        for (int i = 0 ; i < Size ; ++i){
            stringstream ss;
            for (int j = 0 ; j < NodeNameCountVect[i] ; ++j)
                ss << NodeNameList[NodeNameOffsetVect[i] + j];  //  each char into the stringstream
            ss >> FullStrList[i];   //  stringstream into string for each node name
            ss.str(""); //  This and below clear the contents of the stringstream,
            ss.clear(); //  since the >> operator doesn't clear as it extracts
            //cout << FullStrList[i] << endl;   //  for testing
        }
        delete NodeNameList;    //  master is done with full c-array
        bool IsUnique;  //  flag for breaking from for loop
        stringstream ss;    //  used for a full c-array of unique node names
        for (int i = 0 ; i < Size ; ++i){   //  Loop over EVERY name
            IsUnique = true;
            for (int j = 0 ; j < NumUniqueNodes ; ++j)
                if (FullStrList[i].compare(NodeStrList[j]) == 0){   //  check against list of uniques
                    IsUnique = false;
                    break;
                }
            if (IsUnique){
                NodeStrList[NumUniqueNodes] = FullStrList[i];   //  add unique names so others can be checked against them
                ss << NodeStrList[NumUniqueNodes].c_str();  //  build up a string of all unique names back-to-back
                ++NumUniqueNodes;   //  keep a tally of number of unique nodes
            }
        }
        ss >> NodeListStr;  //  make a string of all unique node names
        NodeListCharLen = NodeListStr.size();   //  char length of all unique node names
        NodeListLenVect = new int [NumUniqueNodes]; //  list of unique node name lengths
        /*
         * Because Bcast simply duplicates the buffer of the Bcaster to all cores,
         * the buffer needs to be a char* so that the other cores can have a similar
         * buffer prepared to receive. This wouldn't work if we passed string.c_str()
         * as the buffer, becuase the receiving cores don't have string.c_str() to
         * receive into, and even if they did, c_srt() is a method and can't be used
         * that way.
         */
        NodeNameList = new char [NodeListCharLen];  //  even though c_str is used, allocate necessary memory
        NodeNameList = const_cast<char*>(NodeListStr.c_str());  //  c_str() returns const char*, so need to recast
        for (int i = 0 ; i < NumUniqueNodes ; ++i)  //  fill list of unique node name char lengths
            NodeListLenVect[i] = NodeStrList[i].size();
        /*for (int i = 0 ; i < NumUnique ; ++i)
            cout << UniqueNodeStrList[i] << endl;
        MPI::COMM_WORLD.Abort(1);*/
        //delete NodeStrList;   //  Arrays of string don't need to be deallocated,
        //delete FullStrList;   //  I'm guessing becuase of something weird in the string class.
        delete NodeNameCountVect;
        delete NodeNameOffsetVect;
    }
    /*
     * Now we send the list of node names back to all cores
     * so they can group themselves appropriately.
     */

    //  Bcast the number of nodes in use
    MPI::COMM_WORLD.Bcast(&NumUniqueNodes, 1, MPI::INT, MASTER);
    //  Bcast the full length of all node names
    MPI::COMM_WORLD.Bcast(&NodeListCharLen, 1, MPI::INT, MASTER);

    //  prepare buffers for node name Bcast's
    if (Rank > MASTER){
        NodeListLenVect = new int [NumUniqueNodes];
        NodeNameList = new char [NodeListCharLen];
    }

    //  Lengths of node names for navigating c-string
    MPI::COMM_WORLD.Bcast(NodeListLenVect, NumUniqueNodes, MPI::INT, MASTER);
    //  The actual full list of unique node names
    MPI::COMM_WORLD.Bcast(NodeNameList, NodeListCharLen, MPI::CHAR, MASTER);

    /*
     * Similar to what master did before, each core (incl master)
     * needs to build an actual list of node names as strings so they
     * can compare the c++ way.
     */
    int Offset = 0;
    NodeStrList = new string[NumUniqueNodes];
    for (int i = 0 ; i < NumUniqueNodes ; ++i){
        stringstream ss;
        for (int j = 0 ; j < NodeListLenVect[i] ; ++j)
            ss << NodeNameList[Offset + j];
        ss >> NodeStrList[i];
        ss.str("");
        ss.clear();
        Offset += NodeListLenVect[i];
        //cout << FullStrList[i] << endl;
    }
    //  Now since everyone has the same list, just check your node and find your group.
    int CommGroup = -1;
    for (int i = 0 ; i < NumUniqueNodes ; ++i)
        if (NodeNameStr.compare(NodeStrList[i]) == 0){
            CommGroup = i;
            break;
        }
    if (Rank > MASTER){
        delete NodeListLenVect;
        delete NodeNameList;
    }
    //  In case process fails, error prints and job aborts.
    if (CommGroup < 0){
        cout << "**ERROR** Rank " << Rank << " didn't identify comm group correctly." << endl;
        IsOk = false;
    }

    /*
     * ======================================================================
     * The above method uses c++ strings wherever possible so that things
     * like node name comparisons can be done the c++ way. I'm sure there's
     * a better way to do this because that was way too many lines of code...
     * ======================================================================
     */

    //  Create node communicators
    NodeComm = MPI::COMM_WORLD.Split(CommGroup, 0);
    NodeSize = NodeComm.Get_size();
    NodeRank = NodeComm.Get_rank();

    //  Group for master communicator
    int MasterGroup;
    if (NodeRank == MASTER)
        MasterGroup = 0;
    else
        MasterGroup = MPI_UNDEFINED;

    //  Create master communicator
    MasterComm = MPI::COMM_WORLD.Split(MasterGroup, 0);
    MasterRank = -1;
    MasterSize = -1;
    if (MasterComm != MPI::COMM_NULL){
        MasterRank = MasterComm.Get_rank();
        MasterSize = MasterComm.Get_size();
    }

    MPI::COMM_WORLD.Bcast(&MasterSize, 1, MPI::INT, MASTER);
    NodeComm.Bcast(&MasterRank, 1, MPI::INT, MASTER);

    return IsOk;
}


int main(int argc, char ** argv)
{
    int rank, ranks, size, times, window;
    int args[3];

    /* start up */
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);

    //TEST CODE
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name(processor_name, &len);
    cout<<"Hello World ! from processor "<< processor_name<<" rank "<<rank<<" out of "<<ranks<<" processors\n";
    cout<<"Sum :"<<add(5,6);
    //TEST CODE END

    /* collect hostnames of all the processes */
    gethostname(hostname, sizeof(hostname));
    hostnames = (char*) malloc(sizeof(hostname)*ranks);
    MPI_Gather(hostname, sizeof(hostname), MPI_CHAR, hostnames, sizeof(hostname), MPI_CHAR, 0, MPI_COMM_WORLD);

    /* set job parameters, read values from command line if they're there */
    size = 4096*4;
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
        printf("START mpiGraph v%s\n", VERS);
        printf("MsgSize\t%d\nTimes\t%d\nWindow\t%d\n",size,times,window);
        printf("Procs\t%d\n\n",ranks);
    }

    /* synchronize, then start the run */
    MPI_Barrier(MPI_COMM_WORLD);
    code(rank, ranks, size, times, window);
    MPI_Barrier(MPI_COMM_WORLD);

    /* print memory usage */
/*
  if(rank == 0) { printf("\n"); }
  print_mpi_resources();
*/

    /* mark end of output */
    if (rank == 0) { printf("END mpiGraph\n"); }

    // Test Rank allocation code
    MPI::Intracomm NodeComm;
    MPI::Intracomm MasterComm;

    int NodeRank, MasterRank, NodeSize, MasterSize;
    string NodeNameStr;
    bool b = CommByNode(NodeComm, MasterComm, NodeRank, MasterRank, NodeSize, MasterSize, NodeNameStr);

    cout<<"Bool value:"<<b;
    /* shut down */
    MPI_Finalize();
    return 0;
}