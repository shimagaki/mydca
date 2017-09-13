# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
L, q = 2, 21 
X = np.random.choice(q,L) # Chose L-variables from q-state .
J = np.zeros((L*L, q*q))
J0 = 1.0 

T_equil = 10000
T_relax = 100
N_sampling = 1.0e4 
    
def J_init():
    global J_vec
    global J
    
    fname_in = "J_model_estimator_L"+str(L)+".dat"
    #fname_in = "J_model_estimator_L"+str(L)+"_CD1.dat"
    fin = open(fname_in, "r")
    J_vec = [] 
    i = 0
    for line in fin:
        item = line.split(' ')
        del item[-1] 
        J_vec.append(map(float, item)[0]) 
        i += 1
    ###NOTE############################NOTE
    J = np.reshape( np.copy(J_vec), (L**2,q**2) )
    fin.close()

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
        E_tot = E_local(i)
    return E_tot / 2.0

def Metropolis(i):
    global X 
    E_i = E_local(i)
    Xi_save = X[i]
    X[i] = (X[i]+np.random.choice(q-1)+1) % q  
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
    J_init()
    fname  = "genmodel_data_L"+str(L)+".dat"
    fout = open(fname,"w")
    t = 0 
    while(t<T_equil):
        MonteCarlo_sweep()
        t += 1 
    t = 0 
    while(t<T_relax*N_sampling):
        MonteCarlo_sweep()
        if(t%T_relax == 0):
            for i in range(L):
                fout.write(str(X[i])+" ")
            fout.write("\n")
        t += 1 
    fout.close() 
