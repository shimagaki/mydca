# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
#----------#
q, L = 3, 4 
data_vec = []
model_vec = []

def data_model():
    global data_vec
    global model_vec
    data_vec = [] 
    model_vec = [] 
    
    #fname_data = "Psi_im_data_pdf_CD_L"+str(L)+".dat"
    #fname_model =  "Psi_im_model_pdf_CD_L"+str(L)+".dat"
    #fname_data = "Psi_data_pdf_CD_L"+str(L)+".dat"
    #fname_model =  "Psi_model_pdf_CD_L"+str(L)+".dat"
    
    fname_data = "Correlation_Teacher_L"+str(L)+".dat"
    fname_model =  "Correlation_CD_L"+str(L)+"_hidden.dat"
    
    fin_data = open(fname_data,"r")
    fin_model = open(fname_model,"r")

    i = 0
    for line in fin_data:
        line.replace('\n','')
        data_vec.append(float(line)) 
        i += 1
    i = 0
    for line in fin_model:
        line.replace('\n','')
        #J_model_vec.append(18.0*float(line))
        model_vec.append(float(line))
        i += 1
    fin_data.close()
    fin_model.close()
    
def visualize_estimator():
    n_parameter = len(data_vec)
    x = np.linspace(-0.0,0.06,n_parameter)
    y = np.linspace(-0.0,0.06,n_parameter)
    plt.plot(x,y)
    plt.scatter(data_vec,model_vec)
    #plt.xlabel("true parameter" ,size=18)
    #plt.ylabel("model parameter",size=18)
    #plt.xlabel("$<\delta_{ia}x>_{p(x|a;W)}$",size=28)
    #plt.ylabel(" $<\delta_{ia}x>_{p(x,a|W)}$",size=28)
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    data_model()
    visualize_estimator()

