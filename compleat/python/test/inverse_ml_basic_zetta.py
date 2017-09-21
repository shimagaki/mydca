# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
import time 
q, L = 3, 2 
#----------#
A = []
Psi = np.zeros((L*L,q*q))
Psi_model = np.zeros((L*L,q*q))
#----------#

#----------#
X = np.random.choice(q,L) # Chose L-variables from q-state .
J = np.ones((L*L, q*q)) 
J_data = np.zeros((L*L, q*q))
J_data_vec = [] 
J0 = 1.0 
#----------#
#----------------------------------#
T_equil = 10000
T_relax = 100
M = 1000 # number of sample of average by statistical model. 
eps = 0.05
lr = 0.05 
epoch = 100
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
#----------------------------------#

#----------------------------------#
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
            for j in range(i+1,L):
                a, b = A[i], A[j]
                Psi[i*L+j][a*q+b] += 1.0
                Psi[j*L+i][b*q+a] += 1.0
        l += 1
    Psi = Psi / l
    fin.close()

def average_model():
    global Psi_model
    global X
    Psi_model = np.zeros((L*L,q*q))
    X = np.random.choice(q,L) 
    
    for t in range(T_equil):
        MonteCarlo_sweep()
    acceptances = 0
    for t in range(T_relax*M):
        acceptances += MonteCarlo_sweep()
        if(t%T_relax == 0):
            for i in range(L):
                for j in range(i+1,L):
                    a, b =  X[i],X[j] 
                    Psi_model[i*L+j][a*q+b] += 1.0/M 
                    Psi_model[j*L+i][b*q+a] += 1.0/M
    return acceptances
#----------------------------------#
def gradient_descent():
    error = .0
    global J
    for i in range(L):
        for j in range(i+1,L):
            for a in range(q):
                for b in range(q):
                    dPsi = (Psi[i*L+j][a*q+b] -  Psi_model[i*L+j][a*q+b])
                    error += dPsi**2 
                    J[i*L+j][a*q+b] +=  lr * (dPsi - eps*J[i*L+j][a*q+b]) 
                    
                    dPsi = (Psi[j*L+i][b*q+a] - Psi_model[j*L+i][b*q+a])
                    error += dPsi**2 
                    J[j*L+i][b*q+a] +=  lr * (dPsi - eps*J[j*L+i][b*q+a]) 
    error = np.sqrt(error) / (q**2 * L**2)
    return error
#----------------------------------#
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
    plt.show()

def outputestimator():
    fname = "J_model_estimator_L"+str(L)+".dat"
    f = open(fname, "w")
    J_model_vec = np.reshape( np.copy(J), q**2*L**2)
    for l in range(len(J_model_vec)):
        f.write(str(J_model_vec[l])+" \n")
    f.close()
 
def output_statisticalmodel():
    fname_model = "model_pdf_L"+str(L)+".dat"
    f_model = open(fname_model, "w")
    fname_data = "data_pdf_L"+str(L)+".dat"
    f_data = open(fname_data, "w")
    for i in range(L):
        for j in range(L):
            for a in range(L):
                for b in range(L):
                    f_model.write(str(Psi_model[i*L+j][a*q+b])+"\n")
                    f_data.write(str(Psi[i*L+j][a*q+b])+"\n")
    f_model.close()
    f_data.close()
    
   
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

if __name__ == "__main__":
    #-----------#
    #J_simple_potts()    # initialize of J_data. 
    J_data_init() 
    J_model_init()
    print "J=\n", J
    print "J_data=\n", J_data
    #-----------#
    average()
    #-----------#
    time_start = time.time()
    for t in range(epoch):
        average_model() 
        error = gradient_descent()
        print t, error
    #-----------#
    time_end = time.time()
    dtime = time_end - time_start
    print "#ML: computational time = ", dtime

    outputestimator()
    output_statisticalmodel()
