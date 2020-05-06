//
// Created by prash on 2/21/2020.
//

using namespace std;

#ifndef CS698G_2019_20_2_CLOUDMPI_H
#define CS698G_2019_20_2_CLOUDMPI_H
#endif //CS698G_2019_20_2_CLOUDMPI_H


double ** code(int mypid, int nnodes, long size, long times, long window);

bool CommByNode(MPI::Intracomm &NodeComm, MPI::Intracomm &MasterComm,
                int &NodeRank, int &MasterRank, int &NodeSize, int &MasterSize,
                string &NodeNameStr, MPI::Intracomm &l1_NodeComm);
double ** getDist();

void l2_create_comm(MPI::Intracomm &NodeComm, MPI::Intracomm &MasterComm, MPI_Comm &root_comm, MPI::Intracomm &l1_NodeComm);

bool l1_CommByDatacenter(MPI::Intracomm &NodeComm, MPI::Intracomm &MasterComm,
                   int &NodeRank, int &MasterRank, int &NodeSize, int &MasterSize,
                   string &NodeNameStr, double ** dist, double THRESHOLD);

void l1_create_comm(MPI::Intracomm &NodeComm, MPI::Intracomm &MasterComm, MPI_Comm &root_comm, double ** dist, double THRESHOLD);

// ############################ Custom MPI Collective Calls #####################################
void init();

int MPI_CustomBcast(void *data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator);

int MPI_CustomScatter(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count,
                      MPI_Datatype recv_datatype, int root, MPI_Comm communicator);

int MPI_CustomAllreduce(void* send_data, void* recv_data, int count,
                        MPI_Datatype datatype, MPI_Op op, MPI_Comm communicator);


// ############################ Test (Nothing to look beyond this line) ##########################
int add(int a, int b);