#!/bin/bash
sudo sed -i '/test/d' ./../../etc/hosts
sudo sed -i '1i40.121.90.58 test1useast' /etc/hosts
sudo sed -i '1i52.187.77.210 test2asia' /etc/hosts
sudo sed -i '1i13.76.221.125 test3asia' /etc/hosts
echo "ETC/HOSTS Modified"
rm -rf *
git clone https://github.com/prashplus/CloudMPI
cd CloudMPI
touch hosts
echo "40.121.90.58" >> hosts
echo "52.187.77.210" >> hosts
echo "13.76.221.125" >> hosts
echo "MPI hostfile Modified"
export CC=/usr/bin/mpicc
export CXX=/usr/bin/mpicxx
cmake ../CloudMPI
make 
echo "Make Complete"