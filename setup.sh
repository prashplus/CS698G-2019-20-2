#!/bin/bash
ssh prashant@test1useast.eastus.cloudapp.azure.com 'bash -s' < init.sh
ssh prashant@test2asia.southeastasia.cloudapp.azure.com 'bash -s' < init.sh
ssh prashant@test3asia.southeastasia.cloudapp.azure.com 'bash -s' < init.sh