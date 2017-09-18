# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
import time 
L, q = 2, 3 
#----------#
A = []
Psi = np.zeros((L*L,q*q))
Psi_h = np.zeros((L,q))
Psi_model = np.zeros((L*L,q*q))
Psi_model_h = np.zeros((L,q))

#----------#
X = np.random.choice(q,L) # Chose L-variables from q-state .
J = np.ones((L*L, q*q)) 
J_data = np.zeros((L*L, q*q))
J_data_vec = [] 
J0, h0 = 0.1, 0.1
h = np.ones((L, q)) * h0 
h_data = np.zeros((L, q))
h_data_vec = [] 
J0 = 1.0 

#-----------#
T_equil = int(1e4) 
T_relax = 100
M = 300 # number of sample of average by statistical model. 
eps = 0.05
lr_J, lr_h = 0.5, 0.5 
#----------#
K = 1 # CD_k

def J_data_init():
    global J_data_vec
    global J_data
    
    fname_in = "J_data_L"+str(L)+".dat"
    #fname_in = "test_training_data.dat"
    fin = open(fname_in, "r")
    J_data_vec = [] 
    i = 0
    for line in fin:
        item = line.split(' ')
        del item[-1] 
        J_data_vec.append(map(float, item)[0]) 
        i += 1
    J_data = np.reshape( np.copy(J_data_vec), (L**2,q**2) )
    fin.close()

def h_data_init():
    global h_data_vec
    global  h_data
    
    fname_in = "h_data_L"+str(L)+".dat"
    fin = open(fname_in, "r")
    h_data_vec = [] 
    i = 0
    for line in fin:
        item = float(line) 
        h_data_vec.append(item) 
        i += 1
    h_data = np.reshape( np.copy(h_data_vec), (L,q) )
    fin.close()

def J_simple_potts():
    global J_data
    for i in range(L):
        for j in range(i+1,L):
            for a in range(q):
                J_data[i*L+j][a*q+a] = J0 * (1-2*a) 
                J_data[j*L+i][a*q+a] = J0 * (1-2*a) 

def J_model_init():
    global J
    for i in range(L):
        for a in range(q):
            for b in range(q):
                J[i*L+i][a*q+b] = 0 

def h_model_init():
    global h
    for i in range(L):
        for a in range(q):
            h[i][a] = 0 

def E_local(i):
    a = X[i]
    E_i = .0
    for j in range(L): # NOTE: J[i*L+j][xx] = 0
        b = X[j] 
        E_i += -J[i*L+j][a*q+b] 
    E_i += -h[i][a] 
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
    #X[i] = np.random.choice(q)
    X[i] = (X[i]+np.random.choice(q-1)+1) % q 
    E_i_trial = E_local(i)
    dE = E_i_trial - E_i
    w = np.exp(-dE)
    flag = False
    if(np.random.uniform() < w):
        flag = True
    else:
        X[i] = Xi_save
    return flag 

def MonteCarlo_sweep():
    accepted = .0
    for l in range(L):
        i = np.random.choice(L)
        if(Metropolis(i)):
            accepted +=1
    return accepted

#-----------#
def average():
    global Psi 
    fname_in = "test_training_data_L"+str(L)+".dat"
    fin = open(fname_in, "r")
    l = 0 
    for line in fin:
        item = line.split(' ')
        del item[-1] 
        A = np.copy(map(int, item)) 
        for i in range(L):
            a = A[i]
            for j in range(i+1,L):
                b = A[j]
                Psi[i*L+j][a*q+b] += 1.0
                Psi[j*L+i][b*q+a] += 1.0
        l += 1
    Psi = Psi / l
    fin.close()
    
    global Psi_h 
    fname_in = "test_training_data_L"+str(L)+".dat"
    fin = open(fname_in, "r")
    l = 0 
    for line in fin:
        item = line.split(' ')
        del item[-1] 
        A = np.copy(map(int, item)) 
        for i in range(L):
            a = A[i]
            Psi_h[i][a] += 1.0
        l += 1
    Psi_h = Psi_h / l
    fin.close()

def average_model_CD():
    global Psi_model, Psi_model_h
    global X 
    
    Psi_model = np.zeros((L*L,q*q))
    Psi_model_h = np.zeros((L,q))
    fname_in = "test_training_data_L"+str(L)+".dat"
    fin = open(fname_in, "r")
    
    l = 0 
    for line in fin:
        item = line.split(' ')
        del item[-1]
        X = np.copy(map(int, item))
        
        for k in range(K):
            MonteCarlo_sweep() # n state to n state.
        for i in range(L):
            a =  X[i]
            Psi_model_h[i][a] += 1.0
            for j in range(i+1,L):
                b = X[j] 
                Psi_model[i*L+j][a*q+b] += 1.0
                Psi_model[j*L+i][b*q+a] += 1.0
        l += 1
    Psi_model = Psi_model / l
    Psi_model_h = Psi_model_h / l
    fin.close()

#--------------#
def gradient_descent():
    error, error_h = .0,.0
    global J, h
    for i in range(L):
        for j in range(i+1,L):
            for a in range(q):
                for b in range(q):
                    dPsi = (Psi[i*L+j][a*q+b] -  Psi_model[i*L+j][a*q+b])
                    error += dPsi**2 
                    J[i*L+j][a*q+b] +=  lr_J * (dPsi - eps*J[i*L+j][a*q+b]) 
                    
                    dPsi = (Psi[j*L+i][b*q+a] - Psi_model[j*L+i][b*q+a])
                    error += dPsi**2 
                    J[j*L+i][b*q+a] +=  lr_J * (dPsi - eps*J[j*L+i][b*q+a]) 
    
    for i in range(L):
        for a in range(q):
            dPsi_h = (Psi_h[i][a]-Psi_model_h[i][a])
            error_h += dPsi_h**2
            h[i][a] += lr_h * (dPsi_h - eps*h[i][a])

    error = np.sqrt(error) / (q**2 * L**2)
    error_h = np.sqrt(error_h) / (q*L)
    
    return (error, error_h) 
#----------#
def visualize_estimator():
    n_parameter =  L**2 * q**2
    J_data_vec = np.copy( np.reshape(J_data, n_parameter))    
    J_model_vec = np.copy( np.reshape(J, n_parameter))    
    x = np.linspace(-2,2,n_parameter)
    y = np.linspace(-2,2,n_parameter)
    plt.plot(x,y)
    plt.scatter(J_data_vec,J_model_vec)
    plt.xlabel("true parameter")
    plt.ylabel("model parameter")
    plt.grid(True)
    ##external field is also written isn dehtis wahy # 
    plt.show()

def outputestimator():
    # Coupling Constant 
    fname = "J_model_estimator_CD1_L"+str(L)+".dat"
    f = open(fname, "w")
    J_model_vec = np.reshape( np.copy(J), q**2*L**2)
    for l in range(len(J_model_vec)):
        f.write(str(J_model_vec[l])+" \n")
    f.close()
    
    # Magnetic Field 
    fname = "h_model_estimator_CD1_L"+str(L)+".dat"
    f = open(fname, "w")
    h_model_vec = np.reshape( np.copy(h), q*L)
    for l in range(len(h_model_vec)):
        f.write(str(h_model_vec[l])+" \n")
    f.close()

if __name__ == "__main__":
    #-----------#
    J_data_init()
    h_data_init()
    J_model_init()
    h_model_init()
    print "J=\n", J
    print "J_data=\n", J_data
    print "h=\n", h 
    print "h_data=\n", h_data
    #-----------#
    average()
    #-----------#
    time_start = time.time()
    for t in range(30):
        average_model_CD() 
        error, error_h = gradient_descent()
        print t, error, error_h
    #-----------#
    time_end = time.time()
    dtime = time_end - time_start
    print "#ML: computational time = ", dtime
    outputestimator()
