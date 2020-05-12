#!/bin/bash
zip -r cloudmpi.zip CloudMPI
ssh prashant@test1useast.eastus.cloudapp.azure.com rm -rf *
scp -r cloudmpi.zip prashant@test1useast.eastus.cloudapp.azure.com:cloudmpi.zip
echo "VM 1 ready"
ssh prashant@test2asia.southeastasia.cloudapp.azure.com rm -rf *
scp -r cloudmpi.zip prashant@test2asia.southeastasia.cloudapp.azure.com:cloudmpi.zip
echo "VM 2 ready"
ssh prashant@test3asia.southeastasia.cloudapp.azure.com rm -rf *
scp -r cloudmpi.zip prashant@test3asia.southeastasia.cloudapp.azure.com:cloudmpi.zip
echo "VM 3 ready"
ssh prashant@test4useast.eastus.cloudapp.azure.com rm -rf *
scp -r cloudmpi.zip prashant@test4useast.eastus.cloudapp.azure.com:cloudmpi.zip
echo "VM 4 ready"
bash CloudMPI/setup/setup.sh