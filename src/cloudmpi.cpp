//
// Created by prash on 2/21/2020.
//

#include <bits/stdc++.h>
#include <mpi.h>
#include <sys/resource.h>
#include <unistd.h>

using namespace std;

char  hostname[256];
char* hostnames;
char VERS[] = "1.5";
#define MASTER 0
//#define THRESHOLD 500.0


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

#define __TIME_START__    (g_timeval__start    = MPI_Wtime())
#define __TIME_END_SEND__ (g_timeval__end_send = MPI_Wtime())
#define __TIME_END_RECV__ (g_timeval__end_recv = MPI_Wtime())
#define __TIME_USECS_SEND__ ((g_timeval__end_send - g_timeval__start) * 1000000.0)
#define __TIME_USECS_RECV__ ((g_timeval__end_recv - g_timeval__start) * 1000000.0)
double g_timeval__start, g_timeval__end_send, g_timeval__end_recv;


#endif /* of USE_GETTIMEOFDAY */


MPI::Intracomm l1_NodeComm;
MPI::Intracomm l1_MasterComm;
MPI_Comm l1_root_comm;
MPI::Intracomm l2_NodeComm;
MPI::Intracomm l2_MasterComm;
MPI_Comm l2_root_comm;
double **dist;


/* ****************************************** Network Matrix FUNCTION***********************************************
 * */
double ** code(int mypid, int nnodes, long size, long times, long window)
{
    /* arguments are:
     *   mypid  = rank of this process
     *   nnodes = number of ranks
     *   size   = message size in bytes
     *   times  = number of times to measure bandwidth between task pairs
     *   window = number of outstanding sends and recvs to a single rank
     */
    int i, j, k, w;

    /* collect hostnames of all the processes */
    gethostname(hostname, sizeof(hostname));
    //cout<<hostname;
    hostnames = (char*) malloc(sizeof(hostname)*nnodes);
    MPI_Gather(hostname, sizeof(hostname), MPI_CHAR, hostnames, sizeof(hostname), MPI_CHAR, 0, MPI_COMM_WORLD);
    
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

    if (mypid == 0) {
        printf("%d of %d (100%%)\n", nnodes, nnodes);
        fflush(stdout);
    }
    /* for each node, compute sum of my bandwidths with that node */
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
    double **arr;
    /* rank 0 computes send stats and prints result */
    if (mypid == 0) {
        arr = (double **)malloc(r * sizeof(double *));
        for (i=0; i<r; i++)
            arr[i] = (double *)malloc(r * sizeof(double));
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
//        sendmin /= (double) times;
//        sendsum /= (double) (nnodes)*(nnodes-1)*times;
//        printf("\nSend max\t%f\n", MBsec/sendmin);
//        printf("Send avg\t%f\n", MBsec/sendsum);
//
//        /* print send bandwidth table */
//        printf("\n");
//        printf("Send\t\t\t");
//        for(k=0; k<nnodes; k++) {
//            printf("%s:%d\t", &hostnames[k*sizeof(hostname)], k);
//        }
//        printf("\n");
        for(j=0; j<nnodes; j++) {
            //printf("%s:%d to\t\t", &hostnames[j*sizeof(hostname)], j);
            for(k=0; k<nnodes; k++) {
                double val = allsums[j*nnodes+k];
                arr[j][k]=val/times;
                if (val != 0.0) { val = MBsec * (double) times / val; }
                //printf("%0.3f\t\t", val);
            }
            //printf("\n");
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
//        recvmin /= (double) times;
//        recvsum /= (double) (nnodes)*(nnodes-1)*times;
//        printf("\nRecv max\t%f\n", MBsec/recvmin);
//        printf("Recv avg\t%f\n", MBsec/recvsum);
//
//        /* print receive bandwidth table */
//        printf("\n");
//        printf("Recv\t\t\t");
//        for(k=0; k<nnodes; k++) {
//            printf("%s:%d\t", &hostnames[k*sizeof(hostname)], k);
//        }
//        printf("\n");
        for(j=0; j<nnodes; j++) {
            //printf("%s:%d from\t\t", &hostnames[j*sizeof(hostname)], j);
            for(k=0; k<nnodes; k++) {
                double val = allsums[j*nnodes+k];
                arr[j][k]+=val/times;
                if (val != 0.0) { val = MBsec * (double) times / val; }
                //printf("%0.3f\t\t", val);
            }
            //printf("\n");
        }

    }
//    if(mypid == 0){
//        printf("\n");
//        printf("Combined\t\t\t");
//        for(int k=0; k<nnodes; k++) {
//            printf("%s:%d\t", &hostnames[k*sizeof(hostname)], k);
//        }
//        printf("\n");
//        for(int j=0; j<nnodes; j++) {
//            printf("%s:%d from\t\t", &hostnames[j*sizeof(hostname)], j);
//            for(int k=0; k<nnodes; k++) {
//                printf("%0.3f\t\t", arr[j][k]);
//            }
//            printf("\n");
//        }
//    }

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

    return arr;
}

double ** getDist(){
    int world_rank = MPI::COMM_WORLD.Get_rank();
    int world_size = MPI::COMM_WORLD.Get_size();

    MPI_Barrier(MPI_COMM_WORLD);
    double ** dist = code(world_rank, world_size, 1000, 10, 10);
    MPI_Barrier(MPI_COMM_WORLD);

    return dist;
}

/* ****************************************** Level 1 SETUP ***********************************************
*/

bool l1_CommByDatacenter(MPI::Intracomm &NodeComm, MPI::Intracomm &MasterComm,
                int &NodeRank, int &MasterRank, int &NodeSize, int &MasterSize,
                string &NodeNameStr, double **dist, double THRESHOLD)
{
    bool IsOk = true;

    int Rank = MPI::COMM_WORLD.Get_rank();
    int Size = MPI::COMM_WORLD.Get_size();
    int i,j;
    int CommGroup = -1,*rankmark = (int *) malloc(sizeof(int) * Size);


    if(Rank == 0) {
        int temp = 0;
        vector<pair<int,int>> v;
        for(i=0;i<Size;i++){
            rankmark[i] = -1;
        }

        for (i = 0; i < Size; i++) {
            for (j=0;j< Size;j++){
                if(i!=j && dist[i][j]<THRESHOLD){
                    v.push_back(make_pair(i,j));
                }
            }
        }
        for(i=0;i<v.size();i++){
            //printf("v1 : %d, v2 : %d\n", v[i].first,v[i].second);
            if(rankmark[v[i].first]!= -1){
                if(rankmark[v[i].second] != -1){
                }else{
                    rankmark[v[i].second]= rankmark[v[i].first];
                }
            }else{
                if(rankmark[v[i].second] != -1){
                    rankmark[v[i].first] = rankmark[v[i].second];
                }
                else{
                    rankmark[v[i].first]=temp;
                    rankmark[v[i].second]=temp;
                    temp++;
                }
            }
        }
        for(i=0;i<Size;i++){
            if(rankmark[i]==-1)
                rankmark[i]=temp++;
            //printf("Rankmark [%d] : [%d]\n",i,rankmark[i]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(rankmark,Size,MPI_INT,0,MPI_COMM_WORLD);
    CommGroup = rankmark[Rank];
    //printf("\nRank : %d | L1 CommGroup : %d",Rank, CommGroup);
    //  In case process fails, error prints and job aborts.
    if (CommGroup < 0){
        cout << "**ERROR** Rank " << Rank << " didn't identify comm group correctly." << endl;
        IsOk = false;
    }

    //  Create node communicators
    NodeComm = MPI::COMM_WORLD.Split(CommGroup, 0);
    NodeSize = NodeComm.Get_size();
    NodeRank = NodeComm.Get_rank();
    //cout<<"\nRanking Done L1commby datacenter";
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


void l1_create_comm(MPI::Intracomm &NodeComm, MPI::Intracomm &MasterComm, MPI_Comm &root_comm, double ** dist, double THRESHOLD){

    int world_rank = MPI::COMM_WORLD.Get_rank();
    int world_size = MPI::COMM_WORLD.Get_size();

    int NodeRank, MasterRank, NodeSize, MasterSize;
    string NodeNameStr;
    bool b = l1_CommByDatacenter(NodeComm, MasterComm, NodeRank, MasterRank, NodeSize, MasterSize, NodeNameStr,dist, THRESHOLD);


    //Get the group of processes in MPI_COMM_WORLD
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);

    int rootRing[MasterSize],i,j,k;
    int *rankMatrix,*rankData;

    rankMatrix = (int*)malloc(sizeof(int)*3*world_size);
    rankData = (int*)malloc(sizeof(int)*3);
    rankData[0]=world_rank;
    rankData[1]=NodeRank;
    rankData[2]=MasterRank;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(rankData, 3, MPI_INT, rankMatrix, 3, MPI_INT, MPI_COMM_WORLD);

//    if(rank == 0){
//        printf("\n");
//        printf("\t\tRank Data\n");
//        printf("\n");
//        printf("World Rank\tNode Rank\tMaster Rank\n");
//        for(j=0; j<ranks; j++) {
//            for (k = 0; k < 3; k++) {
//                int val = rankMatrix[j * 3 + k];
//                printf("%d\t\t", val);
//            }
//            printf("\n");
//        }
//    }
    j=0;
    for(i=0;i<world_size;i++){
        if(rankMatrix[i*3+1]==0){
            rootRing[j++]=i;
        }
    }

    // Construct a group containing all of the 0 NodeRanks in world_group
    MPI_Group rootGroup;
    MPI_Group_incl(world_group, MasterSize, rootRing, &rootGroup);

    // Create a new communicator based on the group

    MPI_Comm_create_group(MPI_COMM_WORLD, rootGroup, 0, &root_comm);

    int root_rank = -1, root_size = -1;
    // If this rank isn't in the new communicator, it will be
    // MPI_COMM_NULL. Using MPI_COMM_NULL for MPI_Comm_rank or
    // MPI_Comm_size is erroneous
    if (MPI_COMM_NULL != root_comm) {
        MPI_Comm_rank(root_comm, &root_rank);
        MPI_Comm_size(root_comm, &root_size);
    }

//    cout << '\n' << NodeNameStr;
//    printf("\nReal Rank : %d | NodeRank : %d | L1 MasterRank : %d | L1 Root Rank : %d", world_rank, NodeRank, MasterRank, root_rank);
//    printf("\nReal Size : %d | NodeSize : %d | L1 MasterSize : %d | L1 Root Size : %d\n", world_size, NodeSize, MasterSize, root_size);


    free(rankMatrix);
    free(rankData);

    MPI_Group_free(&world_group);
    MPI_Group_free(&rootGroup);
}

/* ****************************************** Level 2 SETUP ***********************************************
 * */

bool CommByNode(MPI::Intracomm &NodeComm, MPI::Intracomm &MasterComm,
                int &NodeRank, int &MasterRank, int &NodeSize, int &MasterSize,
                string &NodeNameStr, MPI::Intracomm &l1_NodeComm)
{
    bool IsOk = true;
//    For Normal Operation
//    int Rank = MPI::COMM_WORLD.Get_rank();
//    int Size = MPI::COMM_WORLD.Get_size();
    int Size = l1_NodeComm.Get_size();
    int Rank = l1_NodeComm.Get_rank();

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
    //cout<<NodeName;
    //  Prepare a vector for character counts of node names
    if (Rank == MASTER)
        NodeNameCountVect = new int [Size];

    //  Gather node name lengths to master to prepare c-array
    // For Normal Operation
    //MPI::COMM_WORLD.Gather(&NodeNameLen, 1, MPI::INT, NodeNameCountVect, 1, MPI::INT, MASTER);
    MPI_Gather(&NodeNameLen, 1, MPI_INT, NodeNameCountVect, 1, MPI_INT, MASTER, l1_NodeComm);

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
    // Normal Operation
    //MPI::COMM_WORLD.Gatherv(NodeName, NodeNameLen, MPI::CHAR, NodeNameList, NodeNameCountVect, NodeNameOffsetVect, MPI::CHAR, MASTER);
    MPI_Gatherv(NodeName, NodeNameLen, MPI_CHAR, NodeNameList, NodeNameCountVect, NodeNameOffsetVect, MPI_CHAR, MASTER, l1_NodeComm);

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
    // Normal Operation
    //MPI::COMM_WORLD.Bcast(&NumUniqueNodes, 1, MPI::INT, MASTER);
    MPI_Bcast(&NumUniqueNodes, 1, MPI::INT, MASTER, l1_NodeComm);
    //  Bcast the full length of all node names
    // Normal Operation
    //MPI::COMM_WORLD.Bcast(&NodeListCharLen, 1, MPI::INT, MASTER);
    MPI_Bcast(&NodeListCharLen, 1, MPI::INT, MASTER, l1_NodeComm);

    //  prepare buffers for node name Bcast's
    if (Rank > MASTER){
        NodeListLenVect = new int [NumUniqueNodes];
        NodeNameList = new char [NodeListCharLen];
    }

    //  Lengths of node names for navigating c-string
    // Normal Operation
    //MPI::COMM_WORLD.Bcast(NodeListLenVect, NumUniqueNodes, MPI::INT, MASTER);
    MPI_Bcast(NodeListLenVect, NumUniqueNodes, MPI::INT, MASTER, l1_NodeComm);
    //  The actual full list of unique node names
    // Normal Operation
    //MPI::COMM_WORLD.Bcast(NodeNameList, NodeListCharLen, MPI::CHAR, MASTER);
    MPI_Bcast(NodeNameList, NodeListCharLen, MPI::CHAR, MASTER, l1_NodeComm);

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
    // Normal Operation
    //NodeComm = MPI::COMM_WORLD.Split(CommGroup, 0);
    NodeComm = l1_NodeComm.Split(CommGroup, 0);
    NodeSize = NodeComm.Get_size();
    NodeRank = NodeComm.Get_rank();

    //  Group for master communicator
    int MasterGroup;
    if (NodeRank == MASTER)
        MasterGroup = 0;
    else
        MasterGroup = MPI_UNDEFINED;

    //  Create master communicator
    // Normal Operation
    //MasterComm = MPI::COMM_WORLD.Split(MasterGroup, 0);
    MasterComm = l1_NodeComm.Split(MasterGroup, 0);
    MasterRank = -1;
    MasterSize = -1;
    if (MasterComm != MPI::COMM_NULL){
        MasterRank = MasterComm.Get_rank();
        MasterSize = MasterComm.Get_size();
    }

    // Normal Operation
    //MPI::COMM_WORLD.Bcast(&MasterSize, 1, MPI::INT, MASTER);
    MPI_Bcast(&MasterSize, 1, MPI::INT, MASTER, l1_NodeComm);
    NodeComm.Bcast(&MasterRank, 1, MPI::INT, MASTER);

    return IsOk;
}



void l2_create_comm(MPI::Intracomm &NodeComm, MPI::Intracomm &MasterComm, MPI_Comm &root_comm, MPI::Intracomm &l1_NodeComm){
    //Normal Operation
//    int world_rank = MPI::COMM_WORLD.Get_rank();
//    int world_size = MPI::COMM_WORLD.Get_size();
    int world_rank = l1_NodeComm.Get_rank();
    int world_size = l1_NodeComm.Get_size();

    int NodeRank, MasterRank, NodeSize, MasterSize;
    string NodeNameStr;
    bool b = CommByNode(NodeComm, MasterComm, NodeRank, MasterRank, NodeSize, MasterSize, NodeNameStr, l1_NodeComm);

    //Normal Operation
//    MPI_Group world_group;
//    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group world_group;
    MPI_Comm_group(l1_NodeComm, &world_group);

    int rootRing[MasterSize],i,j,k;
    int *rankMatrix,*rankData;

    rankMatrix = (int*)malloc(sizeof(int)*3*world_size);
    rankData = (int*)malloc(sizeof(int)*3);
    rankData[0]=world_rank;
    rankData[1]=NodeRank;
    rankData[2]=MasterRank;

    //Normal Operation
//    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Allgather(rankData, 3, MPI_INT, rankMatrix, 3, MPI_INT, MPI_COMM_WORLD);
    MPI_Barrier(l1_NodeComm);
    MPI_Allgather(rankData, 3, MPI_INT, rankMatrix, 3, MPI_INT, l1_NodeComm);

//    if(rank == 0){
//        printf("\n");
//        printf("\t\tRank Data\n");
//        printf("\n");
//        printf("World Rank\tNode Rank\tMaster Rank\n");
//        for(j=0; j<ranks; j++) {
//            for (k = 0; k < 3; k++) {
//                int val = rankMatrix[j * 3 + k];
//                printf("%d\t\t", val);
//            }
//            printf("\n");
//        }
//    }
    j=0;
    for(i=0;i<world_size;i++){
        if(rankMatrix[i*3+1]==0){
            rootRing[j++]=i;
        }
    }

    // Construct a group containing all of the 0 NodeRanks in world_group
    MPI_Group rootGroup;
    MPI_Group_incl(world_group, MasterSize, rootRing, &rootGroup);

    // Create a new communicator based on the group

    //Normal Operation
//    MPI_Comm_create_group(MPI_COMM_WORLD, rootGroup, 0, &root_comm);
    MPI_Comm_create_group(l1_NodeComm, rootGroup, 0, &root_comm);

    int root_rank = -1, root_size = -1;
    // If this rank isn't in the new communicator, it will be
    // MPI_COMM_NULL. Using MPI_COMM_NULL for MPI_Comm_rank or
    // MPI_Comm_size is erroneous
    if (MPI_COMM_NULL != root_comm) {
        MPI_Comm_rank(root_comm, &root_rank);
        MPI_Comm_size(root_comm, &root_size);
    }

//    cout << '\n' << NodeNameStr;
//    printf("\nReal Rank : %d | NodeRank : %d | L2 MasterRank : %d | L2 Root Rank : %d", world_rank, NodeRank, MasterRank, root_rank);
//    printf("\nReal Size : %d | NodeSize : %d | L2 MasterSize : %d | L2 Root Size : %d\n", world_size, NodeSize, MasterSize, root_size);

    //cout<<"Done L2 Comm";
    free(rankMatrix);
    free(rankData);

    MPI_Group_free(&world_group);
    MPI_Group_free(&rootGroup);
}

// ###################################################################################################################
/* ****************************************** MPI Collective Functions *********************************************** */
// ###################################################################################################################

/* ************************************************* INIT ************************************************************
 * */

void init(){
    dist = getDist();
    // Test Rank allocation code
    l1_create_comm(l1_NodeComm,l1_MasterComm, l1_root_comm, dist, 500);
    l2_create_comm(l2_NodeComm,l2_MasterComm, l2_root_comm, l1_NodeComm);
}

/* ********************************************** MPI Custom Bcast ***************************************************
 * */

int MPI_CustomBcast(void *data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator){

    if(MPI_COMM_NULL != l1_root_comm){
        MPI_Barrier(l1_root_comm);
        MPI_Bcast(data, count, datatype, 0, l1_root_comm);
    }

    if(MPI_COMM_NULL != l2_root_comm){
        MPI_Barrier(l2_root_comm);
        MPI_Bcast(data, count, datatype, 0, l2_root_comm);
    }

    MPI_Barrier(l2_NodeComm);
    MPI_Bcast(data, count, datatype, 0, l2_NodeComm);
    return 1;
}


/* ******************************************MPI Custom Scatter FUNCTION***********************************************
 * */
int MPI_CustomScatter(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count,
                      MPI_Datatype recv_datatype, int root, MPI_Comm communicator){
    int world_rank, world_size;
    MPI_Comm_rank(communicator, &world_rank);
    MPI_Comm_size(communicator, &world_size);

//    MPI_Request request1,request2;
//    MPI_Status status1,status2;

    double *send = (double *)send_data;

    if(world_rank == root){
        for(int i = 0;i< world_size ; i++){
            if(i != root)
                MPI_Send(send + i*send_count, send_count, send_datatype, i, 7, communicator);
        }
    }
    else{
        MPI_Recv(recv_data, recv_count, recv_datatype, world_rank, 7, communicator, MPI_STATUS_IGNORE);
    }
    return 1;
}

/* ******************************************MPI Custom All Reduce FUNCTION***********************************************
 * */

int MPI_CustomAllreduce(void* send_data, void* recv_data, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm communicator){

    MPI_Reduce(send_data, recv_data, count, datatype, op, 0, communicator);

    MPI_Barrier(communicator);

    MPI_CustomBcast(recv_data, count, datatype, 0, communicator);

    return 1;
}

/* ******************************************DEMO FUNCTION***********************************************
 * */
int add(int a, int b) {
    return a+b;
}