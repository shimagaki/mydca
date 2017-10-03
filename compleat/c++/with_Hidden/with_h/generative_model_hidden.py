# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt

q, L = 3, 4 
p = L
X = np.random.choice(q,L) # Chose L-variables from q-state .
Xi0, h0 =  0.3, 0.3 
h = np.zeros((L,q))
Xi = np.ones((L*p,q)) * Xi0

N_sampling = int(1e4) 
T_equil = int(1e3)
t_relax = 30

Xi_data = np.ones((L*p,q)) * Xi0
H = np.zeros(p) # Chose L-variables from q-state .
def Xi_set():
    global Xi_data
    fname_in = "Xi_model_estimator_CD_L"+str(L)+"_Hidden.dat"
    #fname_in = "Xi_model_estimator_CD_L"+str(L)+"_Hidden_c++.dat"
    fin = open(fname_in, "r")
     
    a=0; m=0; i=0; 
    for line in fin:
        line = line.replace('\n','') # do not need to use 
        item = line.split(' ')
        
        Xi_data[i*p+m][a]= float(item[0])
        a += 1
        if(a==q):
            a = 0 
            m += 1
            if(m==p):
                m = 0
                i += 1
    fin.close()

def h_init():
    global  h
    #fname_in = "h_model_estimator_L"+str(L)+".dat"
    fname_in = "h_model_estimator_CD_L"+str(L)+"_Hidden.dat"
    #fname_in = "h_model_estimator_CD_L"+str(L)+"_Hidden_c++.dat"
    fin = open(fname_in, "r")
    a=0; i=0 
    for line in fin:
        line = line.replace('\n','')
        item = line.split(' ')
        h[i][a] = float(item[0]) 
        a += 1
        if(a==q):
            a=0
            i+=1
    fin.close()

def E_local_hidden(i,a):
    global Xi_data, h
    E_i = .0
    for m in range(p): # NOTE: J[i*L+j][xx] = 0
        E_i += - Xi_data[i*p+m][a]*H[m]
    E_i += -h[i][a] 
    return E_i

def sampling_hidden():
    global X, H, Xi_data
    x_0 = np.zeros(p)
    for m in range(p):
        x_0[m] = 0.0
        for i in range(L):
            a = int(X[i] ) 
            x_0[m] += Xi_data[i*p+m][a]
        x_0[m] /= L
    H = np.random.normal(x_0,1.0/L,p)    

def sampling_visible():
    global X
    pi = np.zeros(q)
    X = [] 
    for i in range(L):
        
        E_vec = np.zeros(q)
        for a in range(q):
            E_vec[a] =  E_local_hidden(i,a)
        E_min = np.min(E_vec)
        
        pdis = np.zeros(q) 
        for a in range(q):
            pdis[a] = np.exp(- ( E_local_hidden(i,a) - E_min ))

        pdis_sum = np.sum(pdis)
        pdis = pdis / pdis_sum
        # Does this p give effect on p(= #pattern)?
        X.append( int(np.random.choice(np.arange(q),p=pdis)) ) 

if __name__ == "__main__":
    h_init()
    Xi_set()
    fname_Data  = "test_training_data_L"+str(L)+"_Hidden.dat"
    #fname_Data  = "test_training_data_L"+str(L)+"_Hidden_c++.dat"
    fout_Data = open(fname_Data,"w")
      
    t = 0
    while(t<T_equil):
        sampling_hidden()
        sampling_visible()
        t += 1

    #X = np.random.choice(q,L)
    t =  0
    while(t<N_sampling*t_relax):
        sampling_hidden()
        sampling_visible()
        if(t%t_relax == 0):
            for i in range(L):
                fout_Data.write(str(X[i])+" ")
            fout_Data.write("\n")
        t += 1
    fout_Data.close() 

