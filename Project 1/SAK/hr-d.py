#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

filelist=[]
for i in range(0,20):
    filelist.append("star%s.txt" %i)

for fname in filelist:
    data=np.loadtxt(fname, unpack = True)
    X=data[0]
    Y=data[1]
    plt.plot(X,Y,'.')

plt.gca().invert_xaxis()
plt.yscale('log')
plt.xscale('log')
plt.show()
