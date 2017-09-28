# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
import time 
from numpy import linalg as LA

L,q =  4,3
X_trai = []
X_mode = []
Cov_train = []
Cov_model = []

def calc_cov():
    global X_trai, X_mode 
    global Cov_train, Cov_model
    fname_training = "test_training_data_L"+str(L)+"_c++.dat"
    #fname_model = "genmodel_data_L"+str(L)+".dat"
    fname_model = "genmodel_data_CD_L"+str(L)+".dat"
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
            X_trai = A
        elif(l>0):
            X_trai = np.vstack((X_trai,A))
        l += 1
    #NOTE this is wrong!!, each elements of a Covariant matrix are a matrix not scalar.
    
    Cov_train = np.cov(X_trai.T)
    print "np.shape(Cov_train)=",np.shape(Cov_train) 
    
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
            X_mode = A
        elif(l>0):
            X_mode = np.vstack((X_mode,A))
        l += 1
    Cov_model = np.cov(X_mode.T)

def set_X_cov():
    global X_trai
    global Cov_train, Cov_model
    
    ##-----Training Data-----# 
    fname_training = "test_training_data_L"+str(L)+"_c++.dat"
    f_training = open(fname_training, "r")
    
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
            X_trai = A
        elif(l>0):
            X_trai = np.vstack((X_trai,A))
        l += 1
    
    ##-----Covariance-----# 
    Cov_train = np.zeros((L*q, L*q))
    Cov_model = np.zeros((L*q, L*q))
    
    fname_training = "Correlation_CD_L"+str(L)+"_c++.dat"
    fname_model = "Correlation_Teacher_L"+str(L)+"_c++.dat"
    f_training = open(fname_training, "r")
    f_model = open(fname_model, "r")

    #---- Training Data -----#
    l = 0
    for line in f_training:
        line = line.replace('\n','') # do not need to use 
        item = line.split(' ')
        var = l
        # (L*L,q*q).Mat => (L*q,L*q).Mat
        b = var % q; var -= b; var /= q
        a = var % q; var -= a; var /= q
        j = var % L; var -= j; var /= L
        i = var

        Cov_train[i*q+a][j*q+b] = float(item[0])
        l += 1
    
    #---- Model Data -----#
    l = 0
    for line in f_model:
        line = line.replace('\n','') # do not need to use 
        item = line.split(' ')
        var = l
        # (L*L,q*q).Mat => (L*q,L*q).Mat
        b = var % q; var -= b; var /= q
        a = var % q; var -= a; var /= q
        j = var % L; var -= j; var /= L
        i = var

        Cov_model[i*q+a][j*q+b] = float(item[0])
        l += 1


def plot_cov():
    global Cov_train, Cov_model
    global X_trai
    n_parameter = (L*q)**2
    cov_trainig_vec = np.reshape(Cov_train,n_parameter)
    cov_model_vec = np.reshape(Cov_model,n_parameter)
    x = np.linspace(-0.5,0.5,n_parameter)
    y = np.linspace(-0.5,0.5,n_parameter)
    plt.plot(x,y)
    plt.scatter(cov_trainig_vec,cov_model_vec)
    plt.xlabel("Cov(training)")
    plt.ylabel("Cov(model)")
    plt.grid(True)
    plt.savefig("fig_cd.png")
    #plt.show()

def cast_onto_1st_2nd_eigenspace():
    global Cov_train, Cov_model     
    global X_trai
    W_tra, V_tra = LA.eig(Cov_train)
    W_mod, V_mod = LA.eig(Cov_model)
    V_tra2 = np.copy(V_tra[:2]) 
    V_mod2 = np.copy(V_mod[:2])

    #print "np.shape(np.matrix(X_trai).T)=", np.shape(np.matrix(X_trai).T) 
    print "np.shape(X_trai)=", np.shape(X_trai)
    print "np.shape(V_tra2.T)=", np.shape(V_tra2.T)
    X_trai_cast = np.dot(X_trai, V_tra2.T) 
    X_mode_cast = np.dot(X_trai, V_mod2.T)
     
    plt.subplot(121)
    print "len(X_trai_cast[:,0]) = ", len(X_trai_cast[:,0])
    plt.scatter(X_trai_cast[:,0],X_trai_cast[:,1],s=20,edgecolor='')
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Eigenspace(OriginalData)")
    plt.grid(True)
     
    plt.subplot(122)
    print "len(X_mode_cast[:,0]) = ", len(X_mode_cast[:,0])
    plt.scatter(X_mode_cast[:,0],X_mode_cast[:,1],s=20,edgecolor='')
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Eigenspace(StatisticalModelData)")
    plt.grid(True)
    plt.savefig("cast_sample.png")
    plt.show()

if __name__ == "__main__":
    set_X_cov() 
    #calc_cov()
    cast_onto_1st_2nd_eigenspace() 
    plot_cov()
