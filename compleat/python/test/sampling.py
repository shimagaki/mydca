# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
q, L = 2, 2 
X = np.random.choice(q,L) # Chose L-variables from q-state .
J = np.zeros((L*L, q*q))
J0 = 1.0

T_equil = 100
T_relax = 100
N_sampling = 100
def J_simple_potts():
    global J
    for i in range(L):
        for j in range(i,L):
            for a in range(q):
                J[i*L+j][a*q+a] = J0
                J[j*L+i][a*q+a] = J0
    print "J=\n",J

def E_local(i):
    a = X[i]
    E_i = .0
    for j in range(L): # NOTE: J[i*L+j][xx] = 0
        b = X[j] 
        E_i += -J[i*L+j][a*q+b] 
    return E_i

def E():
    E_tot = .0
    for i in range(L):
        E_tot = 0.5 * E_local(i)
    return E_tot

def Metropolis(i):
    global X 
    E_i = E_local(i)
    Xi_save = X[i]
    X[i] = np.random.choice(q)
    E_i_trial = E_local(i)
    dE = E_i_trial - E_i
    w = np.exp(-dE)
    accepted = False
    if(np.random.uniform() < w):
        accepted = True
    else:
        X[i] = Xi_save
    return accepted

def MonteCarlo_sweep():
    accepted = .0
    for l in range(L):
        i = np.random.choice(L)
        if(Metropolis(i)):
            accepted +=1
    return accepted

if __name__ == "__main__":
    J_simple_potts()
    Engy = .0 
    for t in range(T_equil):
        MonteCarlo_sweep()
    for t in range(T_relax*N_sampling):
        MonteCarlo_sweep()
        if(t%T_relax == 0):
            Engy += E()
    Engy /= N_sampling 
    print "Engy=", Engy


