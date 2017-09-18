# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
#----------#
L, q = 2, 3 
J_data_vec = []
J_model_vec = []

def J_data_J_model():
    global J_data_vec
    global J_model_vec
    J_data_vec = [] 
    J_model_vec = [] 
   
    """
    fname_data = "J_data_L"+str(L)+".dat"
    fname_model = "J_model_estimator_ML_L"+str(L)+".dat"
    #fname_model = "J_model_estimator_CD1_L"+str(L)+".dat"
    """
    
    #Ising Gauge 
    fname_data = "J_data_L"+str(L)+"_IsingGauge.dat"
    fname_model = "J_model_estimator_ML_L"+str(L)+"_IsingGauge.dat"
    #fname_model = "J_model_estimator_CD1_L"+str(L)+"_IsingGauge.dat"
    
    fin_data = open(fname_data,"r")
    fin_model = open(fname_model,"r")

    i = 0
    for line in fin_data:
        item = line.split(' ')
        del item[-1] 
        J_data_vec.append(map(float, item)[0]) 
        i += 1
    fin_data.close()

    i = 0
    for line in fin_model:
        item = line.split(' ')
        del item[-1] 
        J_model_vec.append(map(float, item)[0]) 
        i += 1
    fin_model.close()

def h_data_h_model():
    global h_data_vec
    global h_model_vec
    h_data_vec = [] 
    h_model_vec = [] 
    
    """
    fname_data = "h_data_L"+str(L)+".dat"
    fname_model = "h_model_estimator_ML_L"+str(L)+".dat"
    #fname_model = "h_model_estimator_CD1_L"+str(L)+".dat"
    """
    
    #Ising Gauge 
    fname_data = "h_data_L"+str(L)+"_IsingGauge.dat"
    fname_model = "h_model_estimator_ML_L"+str(L)+"_IsingGauge.dat"
    #fname_model = "h_model_estimator_CD1_L"+str(L)+"_IsingGauge.dat"
    
    fin_data = open(fname_data,"r")
    fin_model = open(fname_model,"r")

    i = 0
    for line in fin_data:
        item = float(line) 
        h_data_vec.append(item) 
        i += 1
    fin_data.close()

    i = 0
    for line in fin_model:
        item = float(line) 
        h_model_vec.append(item) 
        i += 1
    fin_model.close()

def visualize_estimator():
    plt.subplot(121)
    n_parameter_J = len(J_data_vec)
    x_J = np.linspace(-2,2,n_parameter_J)
    y_J = np.linspace(-2,2,n_parameter_J)
    plt.plot(x_J,y_J)
    plt.scatter(J_data_vec,J_model_vec)
    plt.title(r"$J$")
    plt.xlabel("true parameter")
    plt.ylabel("model parameter")
    plt.grid(True)

    plt.subplot(122)
    n_parameter_h = len(h_data_vec)
    x_h = np.linspace(-2,2,n_parameter_h)
    y_h = np.linspace(-2,2,n_parameter_h)
    plt.plot(x_h,y_h)
    plt.scatter(h_data_vec,h_model_vec)
    plt.title(r"$h$")
    plt.xlabel("true parameter")
    plt.ylabel("model parameter")
    plt.grid(True)

    plt.show()

if __name__ == "__main__":
    J_data_J_model()
    h_data_h_model()
    visualize_estimator()

