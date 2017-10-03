# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
#----------#
q, L = 3, 16 
J_data_vec = []
J_model_vec = []

def J_data_J_model():
    global J_data_vec
    global J_model_vec
    J_data_vec = [] 
    J_model_vec = [] 

    #fname_data = "J_data_L"+str(L)+".dat"
    #fname_model = "J_model_estimator_L"+str(L)+".dat"
    
    #fname_data = "J_data_L"+str(L)+"_IsingGauge.dat"
    #fname_model = "J_model_estimator_L"+str(L)+"_IsingGauge.dat"
    
    #fname_data = "data_pdf_L"+str(L)+".dat"
    #fname_model =  "model_pdf_L"+str(L)+".dat"
    
    #fname_data = "fij_L"+str(L)+"_c++.dat"
    #fname_model = "pij_L"+str(L)+"_c++.dat"
    
    fname_data = "Correlation_Teacher_L"+str(L)+"_c++.dat"
    fname_model = "Correlation_CD_L"+str(L)+"_c++.dat"
    
    #fname_data = "fi_CD_L"+str(L)+"_c++.dat"
    #fname_model = "pi_CD_L"+str(L)+"_c++.dat"
    
    fin_data = open(fname_data,"r")
    fin_model = open(fname_model,"r")

    i = 0
    for line in fin_data:
        line.replace('\n','')
        J_data_vec.append(float(line)) 
        i += 1
    i = 0
    for line in fin_model:
        line.replace('\n','')
        #J_model_vec.append(18.0*float(line))
        J_model_vec.append(float(line))
        i += 1
    fin_data.close()
    fin_model.close()
    
def visualize_estimator():
    n_parameter = len(J_data_vec)
    x = np.linspace(-0.03,0.03,n_parameter)
    y = np.linspace(-0.03,0.03,n_parameter)
    plt.plot(x,y)
    plt.scatter(J_data_vec,J_model_vec)
    plt.xlabel("Cov(training)" ,size=18)
    plt.ylabel("Cov(model)",size=18)
    #plt.xlabel("$f_{i}$",size=18)
    #plt.ylabel("$p_{i}$",size=18)
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    J_data_J_model()
    visualize_estimator()


