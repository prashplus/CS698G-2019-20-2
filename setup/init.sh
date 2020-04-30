#!/bin/bash
sudo apt install -y zip unzip mpich
ip1=$(dig +short test1useast.eastus.cloudapp.azure.com)
ip2=$(dig +short test2asia.southeastasia.cloudapp.azure.com)
ip3=$(dig +short test3asia.southeastasia.cloudapp.azure.com)
ip4=$(dig +short test4useast.eastus.cloudapp.azure.com)
sudo sed -i '/test/d' /etc/hosts
sudo sed -i "1i$ip2 test2asia" /etc/hosts
sudo sed -i "1i$ip1 test1useast" /etc/hosts
sudo sed -i "1i$ip3 test3asia" /etc/hosts
sudo sed -i "1i$ip4 test4useast" /etc/hosts
echo "ETC/HOSTS Modified"
# rm -rf *
#git clone https://github.com/prashplus/CloudMPI
unzip cloudmpi.zip
cd CloudMPI
touch hosts
echo "$ip2" >> hosts
echo "$ip1" >> hosts
echo "$ip3" >> hosts
echo "$ip4" >> hosts
echo "MPI hostfile Modified"
export CC=/usr/bin/mpicc
export CXX=/usr/bin/mpicxx
cmake ../CloudMPI
make
echo "Make Complete"