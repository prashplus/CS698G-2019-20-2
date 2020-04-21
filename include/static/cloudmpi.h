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

// ############################ Custom MPI Collective Calls ##############################
int MPI_CustomBcast(void *data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator);

int add(int a, int b);