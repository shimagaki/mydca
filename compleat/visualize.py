# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
#----------#
q, L = 3, 2 
J_data_vec = []
J_model_vec = []

def J_data_J_model():
    global J_data_vec
    global J_model_vec
    J_data_vec = [] 
    J_model_vec = [] 

    fname_data = "J_data_L"+str(L)+".dat"
    fname_model = "J_model_estimator__L"+str(L)+".dat"
    fin_data = open(fname_data,"r")
    fin_model = open(fname_model,"r")

    i = 0
    for line in fin_data:
        item = line.split(' ')
        del item[-1] 
        J_data_vec.append(map(float, item)[0]) 
        i += 1

    i = 0
    for line in fin_model:
        item = line.split(' ')
        del item[-1] 
        J_model_vec.append(map(float, item)[0]) 
        i += 1
    fin_data.close()
    fin_model.close()

def visualize_estimator():
    n_parameter = len(J_data_vec)
    #J_data_vec = np.copy( np.reshape(J_data, n_parameter))    
    #J_model_vec = np.copy( np.reshape(J, n_parameter))    
    x = np.linspace(-2,2,n_parameter)
    y = np.linspace(-2,2,n_parameter)
    plt.plot(x,y)
    plt.scatter(J_data_vec,J_model_vec)
    plt.xlabel("true parameter")
    plt.ylabel("model parameter")
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    J_data_J_model()
    visualize_estimator()











