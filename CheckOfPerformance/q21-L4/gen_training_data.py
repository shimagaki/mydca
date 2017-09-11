# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
q, L = 21, 2 
X = np.random.choice(q,L) # Chose L-variables from q-state .
J = np.zeros((L*L, q*q))
J0 = 1.0 

T_equil = 10000
T_relax = 100
N_sampling = 3*1.0e4

def J_simple_potts():
    global J
    for i in range(L):
        for j in range(i+1,L):
            for a in range(q):
                J[i*L+j][a*q+a] = J0 * (1-2*a) 
                J[j*L+i][a*q+a] = J0 * (1-2*a)    
    print "J=\n",J

def J_generalized_potts():
    global J
    for i in range(L):
        for j in range(i+1,L):
            for a in range(q):
                for b in range(q):
                    r = np.random.normal() 
                    J[i*L+j][a*q+b] = r  
                    J[j*L+i][b*q+a] = r 
                    #J[i*L+j][a*q+a] = J0 * (1-2*a) 
                    #J[j*L+i][a*q+a] = J0 * (1-2*a) 
    fname = "J_data_L"+str(L)+".dat"
    f = open(fname, "w")
    J_vec = np.reshape( np.copy(J), q**2*L**2)
    for l in range(len(J_vec)):
        f.write(str(J_vec[l])+" \n")
    print "J=\n",J
    f.close()

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
    #J_simple_potts()
    J_generalized_potts()
    fname  = "test_training_data_L"+str(L)+".dat"
    fout = open(fname,"w")
    for t in range(T_equil):
        MonteCarlo_sweep()

    t = 0
    while(t<T_relax*N_sampling):
        MonteCarlo_sweep()
        if(t%T_relax == 0):
            for i in range(L):
                fout.write(str(X[i])+" ")
            fout.write("\n")
            #Engy += E()
        t += 1
    fout.close() 
    #Engy /= N_sampling 
    #print "Engy=", Engy


