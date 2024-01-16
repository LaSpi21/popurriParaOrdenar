# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 00:55:51 2020

@author: plasp
"""


from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.linear_model import Perceptron, LogisticRegression
from sklearn.metrics import confusion_matrix
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import  multivariate_normal
X1=multivariate_normal.rvs(size=100,mean=[0,7],cov=[[1,0.0],[0.0, 1]])
X2=multivariate_normal.rvs(size=100,mean=[2,4],cov=[[2,0.0],[0.0, 1]])
plt.scatter(X1[:,0],X1[:,1],color='blue')
plt.scatter(X2[:,0],X2[:,1],color='red')
plt.xlim((-3,10))
plt.ylim((-1,10))
X = np.concatenate((X1,X2))
X = np.c_[X, np.array([1 if i<100 else 2 for i  in range(200)])]


logistic_este=LogisticRegression()
logistic_este.fit(X[:,0:2],X[:,2])

x=np.linspace(-3,10,100)
y=np.linspace(0,10,100)
Xtoplot,Ytoplot=np.meshgrid(x,y)
plt.xlim(-3,10)
plt.ylim(0,10)
Z=logistic_este.predict_proba(np.c_[Xtoplot.ravel(), Ytoplot.ravel()])[:,1].reshape(Xtoplot.shape)
plt.contourf(Xtoplot,Ytoplot,Z,alpha=0.6)
# Z= 
# plt.contourf(Xtoplot,Ytoplot,Z= ,levels=[0,0.5,1],alpha=0.6)

plt.colorbar()
plt.scatter(X1[:,0],X1[:,1],color='blue')
plt.scatter(X2[:,0],X2[:,1],color='red')
plt.legend(loc='upper left',framealpha =0.1)
plt.xlabel(r'PER')
plt.ylabel('USG%')


#%%


#%%


confusio=confusion_matrix(X[:,2],logistic_este.predict(X[:,0:2]))
print(confusio)

#%%


LDA=LinearDiscriminantAnalysis(solver='eigen')
LDA.fit(X[:,0:2],X[:,2])

x=np.linspace(-3,8,100)
y=np.linspace(0,10,100)
Xtoplot,Ytoplot=np.meshgrid(x,y)
plt.xlim(-3,8)
plt.ylim(0,10)
Z=LDA.predict_proba(np.c_[Xtoplot.ravel(), Ytoplot.ravel()])[:,1].reshape(Xtoplot.shape)
plt.contourf(Xtoplot,Ytoplot,Z,levels=[0.0,0.1,0.2,0.3,0.4,0.8,0.9,0.95,0.99,0.999,0.99999,0.99999999999],alpha=0.6)
plt.colorbar()
plt.scatter(X1[:,0],X1[:,1],color='blue')
plt.scatter(X2[:,0],X2[:,1],color='red')

#%%
from sklearn.metrics import confusion_matrix
confusio=confusion_matrix(X[:,2],LDA.predict(X[:,0:2]))
print(confusio)

#%%

LDA.decision_function(X[:,0:2]).shape
plt.hist(LDA.decision_function(X[:100,0:2]),color='orange',histtype='stepfilled',label='No playoff este', alpha = 0.5, bins = 20, cumulative=True)
plt.hist(LDA.decision_function(X[100:201,0:2]),color='green',histtype='stepfilled',label='Playoff este', alpha = 0.5, bins = 20, cumulative=True)
plt.axvline(x=0,color='black',label='Frontera de decision')
plt.legend(loc='upper left',framealpha=0.6)
#%%

LDA.decision_function(X[:,0:2]).shape
plt.hist(LDA.decision_function(X[:100,0:2]),color='orange',histtype='stepfilled',label='No playoff este', alpha = 0.5, bins = 20)
plt.hist(LDA.decision_function(X[100:201,0:2]),color='green',histtype='stepfilled',label='Playoff este', alpha = 0.5, bins = 20)
plt.axvline(x=0,color='black',label='Frontera de decision')
plt.legend(loc='upper left',framealpha=0.6)
#%%

#%%
#Otros datos
from scipy.stats import  multivariate_normal
X1=multivariate_normal.rvs(size=100,mean=[0,7],cov=[[1,0.0],[0.0, 1]])
aux1=multivariate_normal.rvs(size=90,mean=[2,4],cov=[[2,0.0],[0.0, 1]])
aux2=multivariate_normal.rvs(size=10,mean=[9,0],cov=[[0.2,0.0],[0.0, 0.1]])
X2=np.zeros(X1.shape)
X2[:,0]=np.append(aux1[:,0],aux2[:,0])
X2[:,1]=np.append(aux1[:,1],aux2[:,1])
plt.scatter(X1[:,0],X1[:,1],color='blue')
plt.scatter(X2[:,0],X2[:,1],color='red')
plt.xlim((-5,10))
plt.ylim((-1,10))
X = np.concatenate((X1,X2))
X = np.c_[X, np.array([1 if i<100 else 2 for i  in range(200)])]

logistic_este=LogisticRegression()
logistic_este.fit(X[:,0:2],X[:,2])

x=np.linspace(-3,10,100)
y=np.linspace(0,10,100)
Xtoplot,Ytoplot=np.meshgrid(x,y)
plt.xlim(-3,10)
plt.ylim(0,10)
Z=logistic_este.predict_proba(np.c_[Xtoplot.ravel(), Ytoplot.ravel()])[:,1].reshape(Xtoplot.shape)
plt.contourf(Xtoplot,Ytoplot,Z,levels=[0,0.5,1],alpha=0.6)
plt.colorbar()
plt.scatter(X1[:,0],X1[:,1],color='blue')
plt.scatter(X2[:,0],X2[:,1],color='red')
plt.legend(loc='upper left',framealpha =0.1)
plt.xlabel(r'PER')
plt.ylabel('USG%')

#%%

from sklearn.linear_model import Perceptron

X1=multivariate_normal.rvs(size=100,mean=[0,7],cov=[[1,0.0],[0.0, 1]])
X2=multivariate_normal.rvs(size=100,mean=[2,4],cov=[[2,0.0],[0.0, 1]])
plt.scatter(X1[:,0],X1[:,1],color='blue')
plt.scatter(X2[:,0],X2[:,1],color='red')
plt.xlim((-3,10))
plt.ylim((-1,10))
X = np.concatenate((X1,X2))
X = np.c_[X, np.array([1 if i<100 else 2 for i  in range(200)])]
print(X)
perce = Perceptron(max_iter=1000)

# En ese caso, la matriz de diseño es simplemente, x
phi = X[:,0:2].copy()
# print(phi)

perce = perce.fit(phi, X[:,2])

w_perce = (perce.coef_/np.linalg.norm(perce.coef_)).T

print(w_perce)


#%%


def plot_clasi(x, t, ws, labels=[], xp=[-1., 1.], thr=0, spines='zero', equal=True):
    """
    Figura con el resultado del ajuste lineal
    """
    assert len(labels) == len(ws) or len(labels) == 0
    
    if len(labels) == 0:
        labels = np.arange(len(ws)).astype('str')
    
    # Agregemos el vector al plot
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    
    xc1 = x[:, t[0] == 1]
    xc2 = x[:, t[0] == 2]
    
    ax.plot(*xc1, 'ob', mfc='None', label='C1')
    ax.plot(*xc2, 'or', mfc='None', label='C2')

    for i, w in enumerate(ws):
        
        # Ploteo vector de pesos
        x0 = 0.5 * (xp[0] + xp[1])
        ax.quiver(0, thr/w[1], w[0], w[1], color='C{}'.format(i+2), scale=10, label=labels[i],
                 zorder=10)

        # ploteo plano perpendicular
        xp = np.array(xp)
        yp = (thr -w[0]*xp)/w[1] 

        plt.plot(xp, yp, '-', color='C{}'.format(i+2))
        
    # Ploteo línea que une centros de los conjuntos
    mu1 = xc1.mean(axis=1)
    mu2 = xc2.mean(axis=1)
    ax.plot([mu1[0], mu2[0]], [mu1[1], mu2[1]], 'o:k', mfc='None', ms=10)

#     ax.set_xlabel('$x_1$')
#     ax.set_ylabel('$x_2$')
    ax.legend(loc=0, fontsize=16)
    if equal:
        ax.set_aspect('equal')
    
    if spines is not None:
        for a in ['left', 'bottom']:
            ax.spines[a].set_position('zero')
        for a in ['top', 'right']:
            ax.spines[a].set_visible(False)
        
#     ax.set_xlim(-7, 3)
#     ax.set_ylim(-4, 4)
    return
#%%

plot_clasi(X[:,0:2].T, X[:,2:3].T, [w_perce,], xp=[-2, 2])


#%%

plt.scatter(X1[:,0],X1[:,1],color='blue')
plt.scatter(X2[:,0],X2[:,1],color='red')
x=np.linspace(-3,10,100)
y=np.linspace(0,10,100)
Xtoplot,Ytoplot=np.meshgrid(x,y)
plt.xlim(-3,10)
plt.ylim(0,10)



#%%

#Cosas mas reales ponele

a0 = 0.3
a1 = 0.5
a2 = 1
xx = np.linspace(0,2,num=100).reshape(1,-1)
X1 = a0 + a1*xx + np.random.normal(loc = 0, scale = 0.2, size=100)
X2 = a0 + a2*xx + np.random.normal(loc = 0, scale = 0.2, size=100)
X1 = np.vstack((xx,X1)).T
X2 = np.vstack((xx,X2)).T
plt.scatter(X1[:,0],X1[:,1],color='blue')
plt.scatter(X2[:,0],X2[:,1],color='red')
plt.xlim((-0.5,2.5))
plt.ylim((0,3))

#%%
X = np.concatenate((X1,X2))
X = np.c_[X, np.array([1 if i<100 else 2 for i  in range(200)])]


logistic_este=LogisticRegression()
logistic_este.fit(X[:,0:2],X[:,2])

x=np.linspace(0,2,100)
y=np.linspace(0,3,100)
Xtoplot,Ytoplot=np.meshgrid(x,y)
plt.xlim(0,2)
plt.ylim(0,3)
Z=logistic_este.predict_proba(np.c_[Xtoplot.ravel(), Ytoplot.ravel()])[:,1].reshape(Xtoplot.shape)
plt.contourf(Xtoplot,Ytoplot,Z,alpha=0.6)
# Z= 
# plt.contourf(Xtoplot,Ytoplot,Z= ,levels=[0,0.5,1],alpha=0.6)

plt.colorbar()
plt.scatter(X1[:,0],X1[:,1],color='blue')
plt.scatter(X2[:,0],X2[:,1],color='red')
plt.legend(loc='upper left',framealpha =0.1)
plt.xlabel(r'PER')
plt.ylabel('USG%')


#%%


#%%


confusio=confusion_matrix(X[:,2],logistic_este.predict(X[:,0:2]))
print(confusio)

#%%


LDA=LinearDiscriminantAnalysis(solver='eigen')
LDA.fit(X[:,0:2],X[:,2])

x=np.linspace(-3,8,100)
y=np.linspace(0,10,100)
Xtoplot,Ytoplot=np.meshgrid(x,y)
plt.xlim(-3,8)
plt.ylim(0,10)
Z=LDA.predict_proba(np.c_[Xtoplot.ravel(), Ytoplot.ravel()])[:,1].reshape(Xtoplot.shape)
plt.contourf(Xtoplot,Ytoplot,Z,levels=[0.0,0.1,0.2,0.3,0.4,0.8,0.9,0.95,0.99,0.999,0.99999,0.99999999999],alpha=0.6)
plt.colorbar()
plt.scatter(X1[:,0],X1[:,1],color='blue')
plt.scatter(X2[:,0],X2[:,1],color='red')

#%%
from sklearn.metrics import confusion_matrix
confusio=confusion_matrix(X[:,2],LDA.predict(X[:,0:2]))
print(confusio)

#%%


for i in Z:
    print(min(i))