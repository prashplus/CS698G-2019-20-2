#!/bin/bash
sudo sed -i '/asia/d' ./../../etc/hosts
sudo sed -i '/useast/d' ./../../etc/hosts
sudo sed -i '1i40.87.3.199 test1useast' /etc/hosts
sudo sed -i '1i13.67.108.122 test2asia' /etc/hosts
sudo sed -i '1i52.187.35.97 test3asia' /etc/hosts
echo "ETC/HOSTS Modified"
rm -rf CS698G-2019-20-2
git clone https://github.com/prashplus/CS698G-2019-20-2
cd CS698G-2019-20-2
touch hosts
echo "40.87.3.199" >> hosts
echo "13.67.108.122" >> hosts
echo "52.187.35.97" >> hosts
echo "MPI hostfile Modified"
mpicxx -o main.out src.cpp
echo "Program ready to execute"