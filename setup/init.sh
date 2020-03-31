#!/bin/bash
sudo sed -i '/test/d' ./../../etc/hosts
sudo sed -i '1i13.82.103.153 test1useast' /etc/hosts
sudo sed -i '1i168.63.239.148 test2asia' /etc/hosts
sudo sed -i '1i168.63.232.21 test3asia' /etc/hosts
echo "ETC/HOSTS Modified"
rm -rf *
git clone https://github.com/prashplus/CloudMPI
cd CloudMPI
touch hosts
echo "13.82.103.153" >> hosts
echo "168.63.239.148" >> hosts
echo "168.63.232.21" >> hosts
echo "MPI hostfile Modified"
export CC=/usr/bin/mpicc
export CXX=/usr/bin/mpicxx
cmake ../CloudMPI
make 
echo "Make Complete"