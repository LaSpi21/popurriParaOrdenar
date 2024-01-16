# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 21:34:11 2021

@author: plasp
"""

import scipy

# plot "Population" vs "Employed"
from pandas import read_csv
from matplotlib import pyplot
# load the dataset
url = 'https://raw.githubusercontent.com/jbrownlee/Datasets/master/longley.csv'
dataframe = read_csv(url, header=None)
data = dataframe.values
# choose the input and output variables
x, y = data[:, 4], data[:, -1]
# plot input vs output
pyplot.scatter(x, y)
pyplot.show()

#%%


# fit a straight line to the economic data
from numpy import arange
from pandas import read_csv
from scipy.optimize import curve_fit
from matplotlib import pyplot
 
# define the true objective function
def objective(x, a, b):
	return a * x**5 + b
 
# load the dataset
url = 'https://raw.githubusercontent.com/jbrownlee/Datasets/master/longley.csv'
dataframe = read_csv(url, header=None)
data = dataframe.values
# choose the input and output variables
x, y = data[:, 4], data[:, -1]
# curve fit
popt, _ = curve_fit(objective, x, y, bounds = ((0,1),(3,25)))
# summarize the parameter values
a, b = popt
print('y = %.10f * x^5 + %.10f' % (a, b))
# plot input vs output
pyplot.scatter(x, y)
# define a sequence of inputs between the smallest and largest known inputs
x_line = arange(min(x), max(x), 1)
# calculate the output for the range
y_line = objective(x_line, a, b)
# create a line plot for the mapping function
pyplot.plot(x_line, y_line, '--', color='red')
pyplot.show()


#%%


import os
import numpy as np
import glob
from matplotlib import pyplot as plt
from lmfit import minimize, Parameters, report_fit
import numpy as np
import glob
from matplotlib import pyplot as plt
# from lmfit import minimize, Parameters, report_fit
from scipy.optimize import least_squares
import math
# %% Residuals
lista = sorted(glob.glob("*.dat"), key=len)


def abrir(archivo):
    return np.loadtxt(archivo, unpack=True)


q, I = abrir(lista[0])

def residual(params, q, I):
    G = params['G']
    R = params['R']
    d = params['d']
    s = params['s']
    # q_1 = (1/R)*pow(3*d/2,(1/2))

    # D = G*pow(q_1,d)*np.exp(-pow(q_1,2)*pow(R,2)/3)

    # model = (q<=q_1) * (G*np.exp((-pow(q,2)*pow(R,2)/3))) + (q>q_1) * D/pow(q,d)

    # return I - model

    q_1 = (1/R)*math.sqrt(d*(3-s)/2)
    D = G*pow(q_1,d)*np.exp(-pow(q_1,2)*pow(R,2)/(3-s))
    #D = (G/pow(R,d))*math.exp(-d/2*pow(d*(3-s)/2, d/2))
    print(D)
    print(s)
    model = (q<=q_1) * (G/pow(q,s))*np.exp(-pow(q,2)*pow(R,2)/(3-s)) + (q>q_1) * D/pow(q,d+s)
    return I - model

def objective(q, s, R, G, d):
    #q_1 = (1/R)*pow(3*d/2,(1/2))
    #d = d-s
    #D = G*pow(q_1,d)*np.exp(-pow(q_1,2)*pow(R,2)/3)

    #model = (q<=q_1) * (G*np.exp((-pow(q,2)*pow(R,2)/3))) + (q>q_1) * D/pow(q,d)
    #return model
    q_1 = (1/R)*math.sqrt(d*(3-s)/2)
    D = G*pow(q_1,d)*np.exp(-pow(q_1,2)*pow(R,2)/(3-s))
    #D = (G/pow(R,d))*math.exp(-d/2*pow(d*(3-s)/2, d/2))
    print(D)
    print(s)
    model = (q<=q_1) * (G/pow(q,s))*np.exp(-pow(q,2)*pow(R,2)/(3-s)) + (q>q_1) * D/pow(q,d+s)
    return model


params = Parameters()
params.add('G', value = 1,  vary = True, min=0)
params.add('R', value = 10e-10,  vary = True, min=10e-20)
params.add('d', value = 1,  vary = True, min=0)
params.add("s", value = 0,  vary = True, min = 0, max = 2.99)


#Filtro los primeros datos
ar = np.where(q > 0.2)
q_new, I_new = q[ar], I[ar]

out = minimize(residual, params, args = (q_new,I_new))
res = out.params.valuesdict()
print(report_fit(out))
yf = I_new - out.residual

q_line = arange(min(q_new), max(q_new), 0.0001)
# calculate the output for the range
I_model = objective(q_line, res["s"], res["R"], res["G"], res["d"])
plt.plot(q_new,I_new)
plt.plot(q_line, I_model)
# plt.plot(q,yf)
# plt.xlim(0.1,1)
plt.ylim(0.3e6, 4e6)
plt.xlabel('q$(nm^{-1})$')
plt.ylabel('Intensity')

#%%cruve fit

import numpy as np
import glob
from matplotlib import pyplot as plt
# from lmfit import minimize, Parameters, report_fit
from scipy.optimize import least_squares
import scipy
from scipy.optimize import curve_fit
from numpy import arange
from pandas import read_csv
from scipy.optimize import curve_fit
from matplotlib import pyplot
import math

lista = sorted(glob.glob("*.dat"), key=len)


def abrir(archivo):
    return np.loadtxt(archivo, unpack=True)

q, y = abrir(lista[0])

ar = np.where(q > 0.2)
q_new, y_new = q[ar], y[ar]

# define the true objective function
def objective(q, s, R, G, d):
    #general
    #q_1 = (1/R)*(((d-s)*(3-s)/2)**(0.5))
    #D = pow(G,(-pow(q_1,2)*pow(R,2)/(3-s))*pow(q_1,d-s))
    #return (q<=q_1)* pow(G/pow(q,s),(pow(-q,2)*pow(R,2)/(3-s))) + (q>q_1) * D/pow(q,d)
    #general pero s+d es d y d siempre es positivo
    q_1 = (1/R)*math.sqrt(d*(3-s)/2)
    D = G*pow(q_1,d)*np.exp(-pow(q_1,2)*pow(R,2)/(3-s))
    #$D = (G/pow(R,d))*math.exp(-d/2*pow(d*(3-s)/2, d/2))
    return (q<=q_1) * (G/pow(q,s))*np.exp(-pow(q,2)*pow(R,2)/(3-s)) + (q>q_1) * D/pow(q,d+s)
    #q_1 = (pow((1/R*(3*d/2)),(1/2))))
    #D = G**(-)
	#return ((q<=q_1) * (G**((q**2*(R**2)/3))) + (q>q_1) * (G**(-d/2)*pow((3*d/2),(d/2))*(1/R))/pow(q,d)))
    #return (G*q**2 + R*q + d)

popt, _ = curve_fit(objective, q_new, y_new, bounds = ((-np.inf, 10e-100, -np.inf, 10e-10),(2.9,np.inf, np.inf, np.inf)))
# summarize the parameter values
s, R, G, d = popt
#print('y = %.10f * x^5 + %.10f' % (G, R, d))
# plot input vs output
pyplot.scatter(q_new, y_new)
# define a sequence of inputs between the smallest and largest known inputs
q_line = arange(min(q_new), max(q_new), 0.0001)
# calculate the output for the range
y_line = objective(q_line, s, R, G, d)
# create a line plot for the mapping function
pyplot.plot(q_line, y_line, '--', color='red')
pyplot.show()