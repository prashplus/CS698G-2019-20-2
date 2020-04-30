import os
import matplotlib.pyplot as plt
import numpy as np
import statistics

with open('data.txt') as f:
    lines = f.readlines()
    z = [line.split()[0] for line in lines]

rows, cols = (50, 1) 
data_set = [[0 for i in range(cols)] for j in range(rows)]

itr1 = 0
itr2 = 25

for i in range(0,len(z)):
    if i%2 == 0:
        data_set[itr1][0]=z[i]
        itr1 += 1
    else:
        data_set[itr2][0]=z[i]
        itr2 += 1

a = []
k = 0
for i in range(10):
    rows = []
    for j in range(5):
        rows.append(np.asarray(data_set[k]).astype(np.float))
        k+=1
    a.append(rows)

datasize = ['1 KB', '10 KB', '100 KB', '1 MB', '10 MB']

plt.plot(datasize, a[0], linewidth=2, color='black', label='4P Cust_Bcast')
plt.plot(datasize, a[5], linewidth=2,color='silver', label='4P MPI_Bcast')

plt.xlabel('Data Size')
plt.ylabel('Time In Seconds')
plt.title('MPI Custom Bcast vs MPI_Bcast')
plt.legend(loc='best')
plt.show()
plt.clf()

plt.plot(datasize, a[1], linewidth=2,color='navy', label='8P Cust_Bcast')
plt.plot(datasize, a[6], linewidth=2,color='blue', label='8P MPI_Bcast')

plt.xlabel('Data Size')
plt.ylabel('Time In Seconds')
plt.title('MPI Custom Bcast vs MPI_Bcast')
plt.legend(loc='best')
plt.show()
plt.clf()


plt.plot(datasize, a[2], linewidth=2,color='red', label='12P Cust_Bcast')
plt.plot(datasize, a[7], linewidth=2,color='tomato', label='12P MPI_Bcast')

plt.xlabel('Data Size')
plt.ylabel('Time In Seconds')
plt.title('MPI Custom Bcast vs MPI_Bcast')
plt.legend(loc='best')
plt.show()
plt.clf()


plt.plot(datasize, a[3], linewidth=2,color='green', label='16P Cust_Bcast')
plt.plot(datasize, a[8], linewidth=2,color='lime', label='16P MPI_Bcast')

plt.xlabel('Data Size')
plt.ylabel('Time In Seconds')
plt.title('MPI Custom Bcast vs MPI_Bcast')
plt.legend(loc='best')
plt.show()
plt.clf()


plt.plot(datasize, a[4], linewidth=2,color='orange', label='20P Cust_Bcast')
plt.plot(datasize, a[9], linewidth=2,color='wheat', label='20P MPI_Bcast')


plt.xlabel('Data Size')
plt.ylabel('Time In Seconds')
plt.title('MPI Custom Bcast vs MPI_Bcast')
plt.legend(loc='best')
plt.show()
plt.clf()

plt.savefig('image.png')