#!/usr/bin/env python3

import numpy as np
import numpy.random
import matplotlib.pyplot as plt

N = 1000000
N_bins = 100
mean = 0.5
std = 1


x = np.random.normal(mean,std,N)
y = np.zeros(N_bins)
max_x = max(x)+0.0001
min_x = min(x)
#print (max(x),min(x))
bin_size = (max_x-min_x)/N_bins
print (bin_size)

for num in x:
	temp = int((num-min_x)/bin_size)
	y[temp] += 1



#plt.plot(x '.', alpha = 0.2)
plt.bar(np.linspace(min_x,max_x,N_bins),y)
plt.show()
