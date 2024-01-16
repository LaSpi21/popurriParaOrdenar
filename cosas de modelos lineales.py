# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 19:23:26 2020

@author: plasp
"""

# Por si alguien corre en python2
from __future__ import division

# Preparamos todo para correr
import numpy as np
from math import *
from matplotlib import pyplot as plt
%matplotlib inline
from scipy.stats import norm, binom, gamma, poisson, multivariate_normal
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GridSearchCV
import sklearn.preprocessing  as pp
from sklearn.model_selection import cross_val_score
import random
#%%

#Ejercicio 1

a0 = -0.3
a1 = 0.5

X = np.linspace(-1,1,num=100).reshape(1,-1)
Y = a0 + a1*X + np.random.normal(loc = 0, scale = 0.2, size=100)
N = np.vstack((X,Y)).T
np.random.shuffle(N)
print(N)
#%%

# plt.scatter(X,Y)
print(X.shape)
print(Y.shape)
#%%
X_train, X_test, t_train, t_test = train_test_split(N[:,0], N[:,1], test_size=0.2)

print(X.shape)
print(X_train.shape)


#%%
# Array para plotear
xx = np.linspace(-1, 1, 100).reshape([-1, 1])

def plot_data_sine(x, t, ax=None):
    if ax is None:
        ax = plt.gca()
    ax.plot(x, t, 'ob', mfc='None', ms=10)
    ax.plot(xx, a0 + a1*xx, 'g-', lw=2, alpha=0.7, label='Ground Truth')
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.legend(loc=0)
    return

plot_data_sine(X_train, t_train)
#%%


#varios degrees


from sklearn.metrics import mean_squared_error as mse

# Creo el pipeline
pp_sine = Pipeline([('poly', PolynomialFeatures()),
                     ('lr', LinearRegression(fit_intercept=False))])

# Inicializo listas
degrees = [1, 2, 3, 5, 7, 9, 11]
coeffs = []
rmse = []
preds = []

# Itero sobre los grados
for d in degrees:
    # Fijo el grado
    pp_sine.named_steps['poly'].degree = d
    
    # Fiteo (y registro los valores de los parámetros)
    pp_sine.fit(X_train, t_train)
    coeffs.append(pp_sine.named_steps['lr'].coef_)
    
    # Obtengo predicciones
    y_train = pp_sine.predict(X_train)
    
    preds.append(pp_sine.predict(xx))
    rmse.append(np.sqrt(mse(t_train, y_train)))


# Veamos cómo evoluciona la métrica con el número de grados de libertad
plt.figure(figsize=(10,4))
plt.plot(degrees, rmse, 'o-r', mfc='None', ms=10, mew=2, label='Train')
plt.xlabel('Grado de features polinomiales')
plt.ylabel('RMSE')
plt.legend(loc=0)

for i, d in enumerate(degrees):
    print(d, rmse[i])
    
#%%

#Crossval

lr=LinearRegression()
loo_cv=cross_val_score(lr,X.T, Y.T,cv=10,scoring='neg_root_mean_squared_error')
errores = -loo_cv
print(errores.mean(), errores.std())


#%%

def cv_multimodel(grados=range(10), cv=11, plot=True):
        
    rsmes = np.zeros(len(grados))
    std_rsmes = np.zeros(len(grados))
    
    for i, grado in enumerate(grados):
        # Crea un pipeline de sklearn con features de grado "grado"
        modelo = Pipeline([('features', PolynomialFeatures(degree=grado)),#x->[1,x,x^2,..,x^grado]
                           ('regression', LinearRegression(fit_intercept=False))
                          ])
        
        # Hace K-folding
        scores = cross_val_score(modelo, X_test, t_test, cv=cv, scoring='neg_root_mean_squared_error')
            
        # Como se usa un score, hay que pasarlo a Error cambiando de signo.
        rsmes[i] = (-scores).mean()
        std_rsmes[i] = (-scores).std()
        
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        # Sin errores
        #ax.semilogy(grados, rsmes, 'o-', mfc='None')
        
        # Con errores
        ax.errorbar(grados, rsmes, std_rsmes, fmt='-o')
        ax.set_yscale('log')
        
        ax.set_xlabel('Grado')
        ax.set_ylabel('mean RMSE')
                
    return grados, rsmes, std_rsmes

grados, errors, errorst = cv_multimodel()
print("El mejor grado es:{}".format(grados[np.argmin(errors)]))

#%%

# Definir el mejor modelo, elegido a partir de CV
modelo = Pipeline([('features', PolynomialFeatures(degree=3)),
                    ('regression', LinearRegression(fit_intercept=False))
                    ])
# Ajustar con mi conjunto de training.
modelo.fit(X_train, t_train)

# Evaluar en el conjunto de test.
y_test = modelo.predict(X_test)
print('El error del mejor modelo en el conjunto de test es {:.2e}.'.format(np.sqrt(mse(t_test, y_test))))

#%%

#GridSearchCV

modelo = Pipeline([('features', PolynomialFeatures()),
                    ('regression', LinearRegression(fit_intercept=False))
                    ])
parameters={'features__degree':[1, 2, 3, 5, 7, 9, 11]}
grilla=GridSearchCV(modelo,parameters,refit=True)
grilla.fit(X.T, Y.T)
print(grilla.best_params_)
elmejor=grilla.best_estimator_
plt.scatter(X,Y)
  
plt.plot(xx,elmejor.predict(xx))

  
  
#%%


#Mediciones pesadas (con diferentes errores explicitos!)


t_train_scaled=t_train/np.sqrt(np.abs(t_train)+0.001)
Phi=np.hstack((1/np.sqrt(np.abs(t_train)+0.001),X_train/np.sqrt(np.abs(t_train)+0.001)))
# Ahora calculemos el producto de phi por su transpuesta y verifiquemos que la forma es la correcta
pp1 = np.dot(Phi.T, Phi)
# y el producto entre Phi y el vector t
yy = np.dot(Phi.T, t_train)
wml = np.linalg.solve(pp1, yy)
print(wml)

#%%
Phi_preds=np.hstack((xx*0.0+1.0,xx))
preds_manual=np.dot(Phi_preds,wml)
# print(t_train.shape)
plt.errorbar(X_train[:,0], t_train[:,0], np.sqrt(np.abs(t_train)+0.001)[:,0], fmt='ob')#'ob', mfc='None', ms=10)
plt.plot(xx, a0 + a1*xx ,'g-', lw=2, alpha=0.7, label='Ground Truth')
plt.xlabel('x')
plt.ylabel('t')
plt.plot(xx, preds_manual, 'r-', label='Prediction')
plt.ylim(-2.5, 3.3)
plt.legend(loc=0)
plt.title('Grado {}'.format(1))


#%%

# Creo el pipeline
pp_sine = Pipeline([('poly', pp.PolynomialFeatures()),
                     ('lr', LinearRegression(fit_intercept=False))])

# Inicializo listas
degrees = [1, 2, 3, 5, 7, 9, 11]
coeffs = []
rmse = []
rmse_test = []
preds = []
# Itero sobre los grados
for d in degrees:
    # Fijo el grado
    pp_sine.named_steps['poly'].degree = d
    
    # Fiteo (y registro los valores de los parámetros)
    error=np.append([0.00001],np.sqrt(np.abs(t_train)+0.001)[1:,0])
    pp_sine.fit(X_train, t_train,lr__sample_weight=list(map(lambda x: 1/x**2,error)))
    coeffs.append(pp_sine.named_steps['lr'].coef_)
    
    # Obtengo predicciones
    y_train = pp_sine.predict(X_train)
    y_test = pp_sine.predict(X_test)
    preds.append(pp_sine.predict(xx))

    # Calculo el RMSE    
    rmse.append(np.sqrt(mse(t_train, y_train)))
    rmse_test.append(np.sqrt(mse(t_test, y_test)))
plt.figure(figsize=(10,4))
plt.plot(degrees, rmse, 'o-r', mfc='None', ms=10, mew=2, label='Train')
plt.plot(degrees, rmse_test, 'o-b', mfc='None', ms=10, mew=2,label='Test')
plt.xlabel('Grado de features polinomiales')
plt.ylabel('RMSE')
plt.yscale('log')
plt.legend(loc='upper right')
for i, d in enumerate(degrees):
    print(d, rmse[i], rmse_test[i])
    
#%%

ncols = 3
nrows = np.int(np.ceil(len(degrees)/ncols))

fig = plt.figure(figsize=(12, 12))
plt.subplots_adjust(hspace=0.4)

for i, d in enumerate(degrees):
    ax = fig.add_subplot(nrows, ncols, i+1)
    ax.errorbar(X_train[:,0], t_train[:,0], np.append([0.00001],np.sqrt(np.abs(t_train)+0.001)[1:,0]), fmt='ob')#'ob', mfc='None', ms=10)
    ax.plot(xx, a0 + a1*xx,'g-', lw=2, alpha=0.7, label='Ground Truth')
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.plot(xx, preds[i], 'r-', label='Prediction')
    ax.set_ylim(-2.5, 3.3)
    ax.set_title('Grado {}'.format(d))

from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
pp_sine = Pipeline([('poly', PolynomialFeatures()),
                     ('rr', Ridge(fit_intercept=False))])
alphas = np.logspace(-7, -2, 6)
degrees = [2, 3, 5, 7, 9, 11, 13, 15, 17]
param_grid = {'poly__degree': degrees,
             'rr__alpha': alphas}

gg = GridSearchCV(pp_sine, param_grid, 
                  cv=[[np.arange(20), np.arange(20)+20],],
                  scoring='neg_root_mean_squared_error')

gg.fit(np.concatenate([X_train, X_test]), np.concatenate([t_train, t_test]))

