#!/bin/bash
sudo sed -i '/test/d' ./../../etc/hosts
sudo sed -i '1i52.187.32.94 test2asia' /etc/hosts
sudo sed -i '1i40.114.31.169 test1useast' /etc/hosts
sudo sed -i '1i52.187.71.18 test3asia' /etc/hosts
echo "ETC/HOSTS Modified"
rm -rf *
git clone https://github.com/prashplus/CloudMPI
cd CloudMPI
touch hosts
echo "52.187.32.94" >> hosts
echo "40.114.31.169" >> hosts
echo "52.187.71.18" >> hosts
echo "MPI hostfile Modified"
export CC=/usr/bin/mpicc
export CXX=/usr/bin/mpicxx
cmake ../CloudMPI
make 
echo "Make Complete"