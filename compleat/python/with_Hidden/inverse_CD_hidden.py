# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
import time 
q, L = 3,4 
p = L
#----------#
Psi = np.zeros((L*p,q))
Psi_model = np.zeros((L*p,q))
Psi_h = np.zeros((L,q))
Psi_model_h = np.zeros((L,q))

#----------#
X = np.random.choice(q,L) # Chose L-variables from q-state .
Xi0, h0 = 0.1, 0.1
h = np.ones((L, q)) * h0 
Xi = np.ones((L*p,q)) * Xi0

#-----------#
T_equil = int(1e3) 
eps = 1e-6 
lr_J, lr_h = 0.01, 0.01 
#epoch = int(1e2) 
epoch = int(1e2) 
#----------#
K = 1 # CD_k
H = np.zeros(p) # Chose L-variables from q-state .

def E_local_hidden(i,a):
    global Xi, h
    E_i = .0
    for m in range(p): # NOTE: J[i*L+j][xx] = 0
        E_i += - Xi[i*p+m][a]*H[m]
    E_i += -h[i][a] 
    return E_i

def sampling_hidden():
    global X, H, Xi
    x_0 = np.zeros(p)
    for m in range(p):
        x_0[m] = 0.0
        for i in range(L):
            a = int(X[i] ) 
            x_0[m] += Xi[i*p+m][a]
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

# This function is called by the gradient-loop. 
# The function average_data/model_Hidden should be written more efficiently. 
def average_data_Hidden():
    global H, Psi, Psi_h, Xi 
    
    Psi = np.zeros((L*p,q))
    Psi_h = np.zeros((L,q))
    
    fname_in = "test_training_data_L"+str(L)+".dat"
    fin = open(fname_in, "r")
    l = 0 

    for line in fin:
        line = line.replace('\n','')
        item = line.split(' ')
        del item[-1]

        A = np.copy( map(int, item) ) 
        x_0 = np.zeros(p) # = hidden variables 
        for m in range(p):
            x_0[m] = 0.0
            for i in range(L):
                a = int(A[i])
                #NOTE ; Xi is updated in every epoch-loop.
                x_0[m] += Xi[i*p+m][a]
            x_0[m] /= float(L) 
        
        for i in range(L):
            a = int(A[i]) 
            Psi_h[i][a] += 1.0
            for m in range(p):
                Psi[i*p+m][a] += x_0[m]
        l += 1
    Psi = Psi / l
    Psi_h = Psi_h / l
    fin.close()
    
def average_model_CD_Hidden():
    global Psi_model, Psi_model_h
    global X, H #NOTE
    fname_in = "test_training_data_L"+str(L)+".dat"
    fin = open(fname_in, "r")
    
    Psi_model = np.zeros((L*p,q))
    Psi_model_h = np.zeros((L,q))
    
    l = 0 
    for line in fin:
        line = line.replace('\n','')
        item = line.split(' ')
        del item[-1]
        
        X = np.copy(map(int, item))
        
        for k in range(K):
            sampling_hidden() # update of H
            sampling_visible() # update of X
        for i in range(L):
            #for a in range(q):
            a = int(X[i]) 
            Psi_model_h[i][a] += 1.0
            for m in range(p):
                Psi_model[i*p+m][a] += H[m]
     
        l += 1
    Psi_model_h = Psi_model_h / l
    Psi_model  = Psi_model / l
    fin.close()

#--------------------#
def gradient_descent_Hidden():
    error, error_h = .0,.0
    global Xi, h, Psi, Psi_model, Psi_h, Psi_model_h
    for i in range(L):
        for m in range(p):
            for a in range(q):
                dPsi = (Psi[i*p+m][a] -  Psi_model[i*p+m][a])
                error += dPsi**2 
                Xi[i*p+m][a] += lr_J * (dPsi - eps*Xi[i*p+m][a]) 
    for i in range(L):
        for a in range(q):
            dPsi_h = (Psi_h[i][a]-Psi_model_h[i][a])
            error_h += dPsi_h**2
            h[i][a] += lr_h * (dPsi_h - eps*h[i][a])

    error = np.sqrt(error) / (L*p*q)
    error_h = np.sqrt(error_h) / (q*L)
    
    return (error, error_h) 

def outputestimator():
    global Xi, h
    # Coupling Constant 
    fname = "Xi_model_estimator_CD_L"+str(L)+"_Hidden.dat"
    f = open(fname, "w")
    for i in range(L):
        for m in range(p):
            for a in range(q):
                f.write(str(Xi[i*p+m][a])+"\n")
    f.close()
    
    # Magnetic Field 
    fname = "h_model_estimator_CD_L"+str(L)+"_Hidden.dat"
    f = open(fname, "w")
    for i in range(L):
        for a in range(q):
            f.write(str(h[i][a])+"\n")
    f.close()
 
def output_statisticalmodel():
    global Psi_model, Psi
    fname_model = "Psi_model_pdf_CD_L"+str(L)+"_hidden.dat"
    f_model = open(fname_model, "w")
    fname_data = "Psi_data_pdf_CD_L"+str(L)+"_hidden.dat"
    f_data = open(fname_data, "w")
    for i in range(L):
        for m in range(p):
            for a in range(q):
                f_model.write(str(Psi_model[i*p+m][a])+"\n")
                f_data.write(str(Psi[i*p+m][a])+"\n")
    f_model.close()
    f_data.close()
    
if __name__ == "__main__":
    average_data_Hidden()
    time_start = time.time()
    for t in range(epoch):
        average_data_Hidden()
        average_model_CD_Hidden() 
        error, error_h = gradient_descent_Hidden()
        print t, error, error_h
    time_end = time.time()
    dtime = time_end - time_start
    print "#CD: computational time = ", dtime
    outputestimator()
    output_statisticalmodel()
