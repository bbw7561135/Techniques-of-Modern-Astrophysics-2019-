#!/usr/bin/env python3

import numpy as np
import numpy.random
import matplotlib.pyplot as plt

N = 100000000
a = 1.0

x = np.random.random(N)
y = np.random.random(N)

r = (x**2.0 + y**2.0)

print (np.sum(r<1.0)*4.0/N)
