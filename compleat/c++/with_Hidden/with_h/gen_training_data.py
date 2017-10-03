# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt

q, L = 3,  4 
X = np.random.choice(q,L) # Chose L-variables from q-state .
J = np.zeros((L*L, q*q))
h = np.zeros((L,q))
J0, h0 =  0.3, 0.3 

T_equil = int(1e3) 
T_relax = int(1e2) 
N_sampling = int(1e4) 

def h_generalized_potts():
    global h
    for i in range(L):
        for a in range(q):
            r = np.random.normal() * h0
            h[i][a] = r 
    fname = "h_data_L"+str(L)+".dat"
    f = open(fname, "w")
    h_vec = np.reshape(np.copy(h),q*L)
    for l in range(len(h_vec)):
        f.write(str(h_vec[l])+ "\n")
    print "h=\n", h
    f.close()

def J_generalized_potts():
    global J
    for i in range(L):
        for j in range(i+1,L):
            for a in range(q):
                for b in range(q):
                    r = np.random.normal() * J0 
                    J[i*L+j][a*q+b] = r  
                    J[j*L+i][b*q+a] = r 
                    #J[i*L+j][a*q+a] = J0 * (1-2*a) 
                    #J[j*L+i][a*q+a] = J0 * (1-2*a) 
    fname = "J_data_L"+str(L)+".dat"
    f = open(fname, "w")
    J_vec = np.reshape( np.copy(J), q**2*L**2)
    for l in range(len(J_vec)):
        f.write(str(J_vec[l])+"\n")
    print "J=\n",J
    f.close()

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
    J_generalized_potts()
    h_generalized_potts()
    fname_J  = "test_training_data_L"+str(L)+".dat"
    fout_J = open(fname_J,"w")
    
    t =  0
    while(t<T_equil):
        MonteCarlo_sweep()
        t += 1
    
    t =  0
    while(t<T_relax*N_sampling):
        MonteCarlo_sweep()
        if(t%T_relax == 0):
            for i in range(L):
                fout_J.write(str(X[i])+" ")
            fout_J.write("\n")
        t += 1
    fout_J.close() 

