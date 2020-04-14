#!/bin/bash
sudo sed -i '/test/d' ./../../etc/hosts
sudo sed -i '1i13.76.32.27 test2asia' /etc/hosts
sudo sed -i '1i52.226.78.96 test1useast' /etc/hosts
sudo sed -i '1i13.76.39.242 test3asia' /etc/hosts
sudo sed -i '1i13.92.101.49 test4useast' /etc/hosts
echo "ETC/HOSTS Modified"
rm -rf *
git clone https://github.com/prashplus/CloudMPI
cd CloudMPI
touch hosts
echo "13.76.32.27" >> hosts
echo "52.226.78.96" >> hosts
echo "13.76.39.242" >> hosts
echo "13.92.101.49" >> hosts
echo "MPI hostfile Modified"
export CC=/usr/bin/mpicc
export CXX=/usr/bin/mpicxx
cmake ../CloudMPI
make 
echo "Make Complete"