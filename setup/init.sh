#!/bin/bash
sudo sed -i '/test/d' ./../../etc/hosts
sudo sed -i '1i23.100.29.83 test1useast' /etc/hosts
sudo sed -i '1i168.63.251.205 test2asia' /etc/hosts
sudo sed -i '1i137.116.151.143 test3asia' /etc/hosts
echo "ETC/HOSTS Modified"
rm -rf *
git clone https://github.com/prashplus/CloudMPI
cd CloudMPI
touch hosts
echo "23.100.29.83" >> hosts
echo "168.63.251.205" >> hosts
echo "137.116.151.143" >> hosts
echo "MPI hostfile Modified"
export CC=/usr/bin/mpicc
export CXX=/usr/bin/mpicxx
cmake ../CloudMPI
make 
echo "Make Complete"