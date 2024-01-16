# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 15:01:13 2020

@author: plasp
"""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import scipy.cluster.vq as scv


def red_black_green():
    cdict = {
       'red': ((0.0, 0.0, 0.0),
               (0.5, 0.0, 0.0),
               (1.0, 1.0, 1.0)),
       'blue': ((0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0)),
       'green': ((0.0, 0.0, 1.0),
                 (0.5, 0.0, 0.0),
                 (1.0, 0.0, 0.0))
       }

    my_cmap = mpl.colors.LinearSegmentedColormap(
        'my_colormap', cdict, 100)

    return my_cmap

A = red_black_green()

def colormap2arr(arr,cmap):    
    # http://stackoverflow.com/questions/3720840/how-to-reverse-color-map-image-to-scalar-values/3722674#3722674
    print(1)
    gradient=cmap(np.linspace(0.0,1.0,100))
    print(2)
    # Reshape arr to something like (240*240, 4), all the 4-tuples in a long list...
    arr2=arr.reshape((arr.shape[0]*arr.shape[1],arr.shape[2]))
    print(3)
    print(arr2.shape)
    print(gradient.shape)
    # Use vector quantization to shift the values in arr2 to the nearest point in
    # the code book (gradient).
    code,dist=scv.vq(arr2,gradient)
    print(4)
    # code is an array of length arr2 (240*240), holding the code book index for
    # each observation. (arr2 are the "observations".)
    # Scale the values so they are from 0 to 1.
    values=code.astype('float')/gradient.shape[0]
    # Reshape values back to (240,240)
    values=values.reshape(arr.shape[0],arr.shape[1])
    values=values[::-1]
    return values

arr=plt.imread('mri_demo3.png')
values=colormap2arr(arr,A)    
# Proof that it works:
plt.imshow(values,interpolation='bilinear', cmap=A,
           origin='lower', extent=[-3,3,-3,3])
plt.show()

def trans(values):
    return (values[:,2]-0.5)*2/.47