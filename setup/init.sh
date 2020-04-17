#!/bin/bash
sudo apt install -y zip unzip
sudo sed -i '/test/d' ./../../etc/hosts
sudo sed -i '1i137.116.139.58 test2asia' /etc/hosts
sudo sed -i '1i13.82.82.247 test1useast' /etc/hosts
sudo sed -i '1i137.116.139.117 test3asia' /etc/hosts
sudo sed -i '1i13.82.83.246 test4useast' /etc/hosts
echo "ETC/HOSTS Modified"
# rm -rf *
#git clone https://github.com/prashplus/CloudMPI
unzip cloudmpi.zip
cd CloudMPI
touch hosts
echo "137.116.139.58" >> hosts
echo "13.82.82.247" >> hosts
echo "137.116.139.117" >> hosts
echo "13.82.83.246" >> hosts
echo "MPI hostfile Modified"
export CC=/usr/bin/mpicc
export CXX=/usr/bin/mpicxx
cmake ../CloudMPI
make 
echo "Make Complete"