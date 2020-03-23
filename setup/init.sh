#!/bin/bash
sudo sed -i '/asia/d' ./../../etc/hosts
sudo sed -i '/useast/d' ./../../etc/hosts
sudo sed -i '1i40.121.156.19 test1useast' /etc/hosts
sudo sed -i '1i52.230.7.177 test2asia' /etc/hosts
sudo sed -i '1i52.230.4.24 test3asia' /etc/hosts
echo "ETC/HOSTS Modified"
rm -rf CloudMPI
git clone https://github.com/prashplus/CloudMPI
cd CloudMPI
touch hosts
echo "40.121.156.19" >> hosts
echo "52.230.7.177" >> hosts
echo "52.230.4.24" >> hosts
echo "MPI hostfile Modified"
cd cmake-build-release
make clean
make
echo "Program ready to execute"