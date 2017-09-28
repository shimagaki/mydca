# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
L, q = 4, 3 
X = np.random.choice(q,L) # Chose L-variables from q-state .
J = np.zeros((L**2, q**2))
h = np.zeros((L,q))
J0 = 1.0 

T_equil = int(1e3) 
T_relax = int(1e2) 
N_sampling = int(1e3) 
    
def J_init():
    global J_vec
    global J
    
    #fname_in = "J_model_estimator_L"+str(L)+".dat"
    fname_in = "J_model_estimator_CD_L"+str(L)+".dat"
    fin = open(fname_in, "r")
    J_vec = [] 
    i = 0
    for line in fin:
        line = line.replace('\n','')
        item = line.split(' ')
        J_vec.append(float(item[0])) 
        i += 1
    ###NOTE############################NOTE
    J = np.reshape( np.copy(J_vec), (L**2,q**2) )
    fin.close()

def h_init():
    global h_data_vec
    global  h
    
    #fname_in = "h_model_estimator_L"+str(L)+".dat"
    fname_in = "h_model_estimator_CD_L"+str(L)+".dat"
    fin = open(fname_in, "r")
    h_vec = [] 
    i = 0
    for line in fin:
        line = line.replace('\n','')
        item = line.split(' ')
        h_vec.append(float(item[0])) 
        i += 1
    h = np.reshape( np.copy(h_vec), (L,q) )
    fin.close()

def E_local(i):
    a = X[i]
    E_i = .0
    for j in range(L): # NOTE: J[i*L+j][xx] = 0
        b = X[j] 
        E_i += -J[i*L+j][a*q+b] 
    E_i += - h[i][a] 
    return E_i

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
    h_init()
    #fname  = "genmodel_data_L"+str(L)+".dat"
    fname  = "genmodel_data_CD_L"+str(L)+".dat"
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
