//
// Created by prash on 2/21/2020.
//

using namespace std;

#ifndef CS698G_2019_20_2_CLOUDMPI_H
#define CS698G_2019_20_2_CLOUDMPI_H

#endif //CS698G_2019_20_2_CLOUDMPI_H

double ** code(int mypid, int nnodes, int size, int times, int window);

bool CommByNode(MPI::Intracomm &NodeComm, MPI::Intracomm &MasterComm,
                int &NodeRank, int &MasterRank, int &NodeSize, int &MasterSize,
                string &NodeNameStr);

void l2_create_comm(MPI::Intracomm &NodeComm, MPI::Intracomm &MasterComm, MPI_Comm &root_comm);

bool l1_CommByDatacenter(MPI::Intracomm &NodeComm, MPI::Intracomm &MasterComm,
                   int &NodeRank, int &MasterRank, int &NodeSize, int &MasterSize,
                   string &NodeNameStr, double ** dist);
void l1_create_comm(MPI::Intracomm &NodeComm, MPI::Intracomm &MasterComm, MPI_Comm &root_comm, double ** dist);

int add(int a, int b);