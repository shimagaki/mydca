# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
#-------#
L, q = 8, 3
var1 = []
var2  = []

def set_variables():
    global var1, var2 
    var1 = [] 
    var2 = [] 
    #------Parameter-------#
    #fname1 = "J_data_L"+str(L)+"_c++.dat"
    #fname2 = "J_model_estimator_L"+str(L)+"_c++.dat"
    #fname2 = "J_model_estimator_L"+str(L)+"_CD1_c++.dat"
    #------Ising-gauge------#
    fname1 = "J_data_L"+str(L)+"_IsingGauge_c++.dat"
    fname2 = "J_model_estimator_L"+str(L)+"_IsingGauge_c++.dat"
    
    #------Covariance-------#
    #fname1 = "Correlation_Teacher_L"+str(L)+"_c++.dat"    
    #fname2 = "Correlation_ML_L"+str(L)+"_c++.dat"    
    
    fin1 = open(fname1,"r")
    fin2 = open(fname2,"r")

    i = 0
    for line in fin1:
        var1.append(float(line)) 
        i += 1
    fin1.close()

    i = 0
    for line in fin2:
        var2.append(float(line)) 
        i += 1
    fin2.close()

def visualize_estimator():
    n_parameter = len(var1)
    x = np.linspace(-1,1,n_parameter)
    y = np.linspace(-1,1,n_parameter)
    plt.plot(x,y)
    plt.scatter(var1,var2)
    plt.xlabel("true parameter")
    plt.ylabel("model parameter")
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    set_variables() 
    visualize_estimator()

