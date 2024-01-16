# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 02:50:52 2020

@author: plasp
"""

from numba import jit
import numpy as np
from timeit import default_timer as timer


x = np.arange(1000000000).reshape(1000, 1000,1000)

@jit(nopython=True, parallel = True, fastmath = True) # Set "nopython" mode for best performance, equivalent to @njit
def go_fast(a): # Function is compiled to machine code when called the first time
    trace = 0.0
    for i in range(a.shape[0]):   # Numba likes loops
        trace += np.tanh(a[i, i, i]) # Numba likes NumPy functions
    return a + trace              # Numba likes NumPy broadcasting

start = timer()
go_fast(x)
print(timer() - start)

start = timer()
go_fast(x)
print(timer() - start)

#%%

from numba import jit
import numpy as np
from timeit import default_timer as timer


x = np.arange(1000000000).reshape(1000, 1000,1000)

def go_fast(a): # Function is compiled to machine code when called the first time
    trace = 0.0
    for i in range(a.shape[0]):   # Numba likes loops
        trace += np.tanh(a[i, i, i]) # Numba likes NumPy functions
    return a + trace              # Numba likes NumPy broadcasting

start = timer()
go_fast(x)
print(timer() - start)

start = timer()
go_fast(x)
print(timer() - start)


    
