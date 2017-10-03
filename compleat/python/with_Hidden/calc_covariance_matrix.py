# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
import time 
L,q =  4,3 
Cov_training = []
Cov_model = []

def calc_cov():
    global Cov_training
    global Cov_model
    fname_training = "test_training_data_L"+str(L)+".dat"
    #fname_model = "genmodel_data_L"+str(L)+".dat"
    #fname_model = "genmodel_data_CD_L"+str(L)+".dat"
    fname_model = "test_training_data_L"+str(L)+"_Hidden.dat"
    f_training = open(fname_training, "r")
    f_model = open(fname_model, "r")
    l = 0 
    for line in f_training:
        line = line.replace('\n','') # do not need to use 
        item = line.split(' ')
        del item[-1]
        B = np.copy(map(int, item))
        #A should be a vector, which has length L*q
        A = np.zeros(L*q) 
        for i in range(L):
            r = B[i]
            s = q*i + r
            A[s] += 1

        if(l==0):
            Data = A
        elif(l>0):
            Data = np.vstack((Data,A))
        l += 1
    #NOTE this is wrong!!, each elements of a Covariant matrix are a matrix not scalar.
    
    Cov_training = np.cov(Data.T)
    
    l = 0 
    for line in f_model:
        line = line.replace('\n','')
        item = line.split(' ')
        del item[-1]
        B = np.copy(map(int, item))
        #A should be a vector, which has length L*q
        A = np.zeros(L*q) 
        for i in range(L):
            r = B[i]
            s = q*i + r
            A[s] += 1

        if(l==0):
            Data = A
        elif(l>0):
            Data = np.vstack((Data,A))
        l += 1
    Cov_model = np.cov(Data.T)

def analysis_cov():
    n_parameter = (L*q)**2
    cov_trainig_vec = np.reshape(Cov_training,n_parameter)
    cov_model_vec = np.reshape(Cov_model,n_parameter)
    fname = "Cov_data_modle_L"+str(L)+".dat"
    fout = open(fname,"w")
    for l in range(n_parameter):
        fout.write(str(cov_trainig_vec[l]) +" "+str(cov_model_vec[l]) + "\n")
    fout.close() 
    
    x = np.linspace(-0.5,0.5,n_parameter)
    y = np.linspace(-0.5,0.5,n_parameter)
    plt.plot(x,y)
    plt.scatter(cov_trainig_vec,cov_model_vec)
    plt.xlabel("Cov(training)")
    plt.ylabel("Cov(model)")
    plt.grid(True)
    plt.savefig("fig_cd.png")
    plt.show()

if __name__ == "__main__":
    calc_cov()
    analysis_cov()

