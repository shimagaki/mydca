# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
import time 

q, L = 3, 4
p = L
N = 100
x=np.linspace(0.0,3.0,N)
y=[]
pdis = [0.1, 0.05, 0.05, 0.2, 0.4, 0.2]
xdis = np.arange(1,len(pdis)+1) 
print xdis 
def accume_y():
    global y
    y0 = 0
    for i in range(N):
        y0 += x[i]
        y.append(y0)
    y_sum = np.sum(y)
    y /= y0 

def originaldist():
    fdis = np.zeros(len(pdis)) 
    for i in range(N):
        a = np.random.choice(xdis,p=pdis) 
        fdis[a-1] += 1 
    fdis /= float(N)
    return fdis

def plot():
    #plt.plot(x,y)
    fdis = originaldist()
    plt.plot(map(float,xdis),fdis, c='red',label="resample")
    plt.plot(map(float,xdis),pdis, c='blue',label="original")
    plt.legend()
    plt.show()




if __name__ == "__main__":
    #accume_y()
    plot()
