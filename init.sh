#!/bin/bash
sudo sed -i '/asia/d' ./../../etc/hosts
sudo sed -i '/useast/d' ./../../etc/hosts
sudo sed -i '1i40.117.169.110 test1useast' /etc/hosts
sudo sed -i '1i52.230.120.181 test2asia' /etc/hosts
sudo sed -i '1i52.187.147.234 test3asia' /etc/hosts
echo "ETC/HOSTS Modified"
rm hosts
touch hosts
echo "40.117.169.110" >> hosts
echo "52.230.120.181" >> hosts
echo "52.187.147.234" >> hosts
echo "MPI hostfile Modified"
rm -rf CS698G-2019-20-2
git clone https://github.com/prashplus/CS698G-2019-20-2
cd CS698G-2019-20-2
mpicc -o main.out src.c
echo "Program ready to execute"