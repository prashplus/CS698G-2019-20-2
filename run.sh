#!/bin/bash
mpiexec -np 4 -f hosts ./main 10000
mpiexec -np 8 -f hosts ./main 10000
mpiexec -np 12 -f hosts ./main 10000
mpiexec -np 16 -f hosts ./main 10000
mpiexec -np 20 -f hosts ./main 10000