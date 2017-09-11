# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
q, L = 3, 3 
A = []
Psi = np.zeros((L*L,q*q))

def read():
    global A
    fname_in = "test_data4.dat"
    fin = open(fname_in, "r")
    i = 0 
    for line in fin:
        item = line.split(' ')
        #del item[-1] 
     
        a = np.copy(map(int, item)) 
        if(i==0):
            A = a 
        if(i==1):
            A = np.append([np.copy(A)],[a],axis = 0)
        if(i>1):
            A = np.append(np.copy(A),[a],axis = 0)
        i += 1
    fin.close()

def average():
    global Psi
    M = len(A)
    for m in range(M):
        for i in range(L):
            for j in range(L):
                if(i!=j):
                    a, b = A[m][i], A[m][j] 
                    Psi[i*L+j][a*q+b] += 1.0/M 
if __name__ == "__main__":
    read()
    average()
    for i in range(L):
        for j in range(L):
            print "\n",i,",",j  
            for a in range(q):
                print Psi[i*L+j][a*q], " ",Psi[i*L+j][a*q+1], " ",Psi[i*L+j][a*q+2]
