# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
L, q = 4, 3 
Jmean = []
Jmean_k = []

def J_init(fJname):
    
    fin = open(fJname, "r")
    J_vec = [] 
    i = 0
    for line in fin:
        line = line.replace('\n','') 
        print "line=", line
        J_vec.append(float(line)) 
        i += 1
    ###NOTE#######################NOTE###
    print "i=", i
    J = np.reshape( np.copy(J_vec), (L**2,q**2) )
    fin.close()
    return J

def h_init(fhname):
    fin = open(fhname, "r")
    h_vec = [] 
    i = 0
    for line in fin:
        line = line.replace('\n','')
        print "lineh=", line
        h_vec.append(float(line)) 
        i += 1
    h = np.reshape( np.copy(h_vec), (L,q) )
    fin.close()
    return h

def Ising_Gauge_J(J):
    global Jmean, Jmean_k
    J_IG = np.copy(J)
    
    Jmean = np.zeros(L**2)
    Jmean_k = np.zeros((L**2,q))
    for i in range(L):
        for j in range(1+i,L):
            Jij_b = np.zeros(q) # sum of col
            Jija_ = np.zeros(q) # sum of row 
                 
            for c in range(q):
                for d in range(q):
                    Jij_b[c] += J[i*L+j][d*q+c]/q
                    Jija_[c] += J[i*L+j][c*q+d]/q
                    Jmean[i*L+j] += J[i*L+j][d*q+c]/q**2
                    Jmean[j*L+i] += J[j*L+i][c*q+d]/q**2  
                Jmean_k[i*L+j][c] =  Jija_[c]
                Jmean_k[j*L+i][c] =  Jij_b[c]
            Jij__ = sum(Jij_b) / q 
            
            for a in range(q):
                for b in range(q):
                    J_IG[i*L+j][a*q+b] = J[i*L+j][a*q+b] \
                            - Jij_b[b] - Jija_[a] + Jij__
    
                    J_IG[j*L+i][b*q+a] = J_IG[i*L+j][a*q+b]
    return J_IG

def Ising_Gauge_h(h):
    h_IG = np.copy(h)
    CorrectJ = np.zeros((L,q)) 
    """
    for i in range(L):
        for j in range(L):
            if(j != i):
                for a in range(q):
                    for b in range(q):
                        CorrectJ[i*L+j][a] += J[i*L+j][a*L+b] / q 
                ave_CorrectJ = sum(CorrectJ[i*L+j,:]) / q
                CorrectJ[i*L+j,:] = np.copy(CorrectJ[i*L+j,:]) - ave_CorrectJ
    
    Correcth = np.zeros((L,q))
    for i in range(L): 
        for j in range(L):
            if(j != i):
                Correcth[i,:] = np.copy(Correcth[i,:]) + CorrectJ[i*L+j,:] 
    """
    for i in range(L):
        for j in range(L):
            if(j!=i):
                for k in range(q):
                    CorrectJ[i][k] +=Jmean_k[i*L+j][k]-Jmean[i*L+j] 
                    CorrectJ[i][k] +=Jmean_k[i*L+j][k]-Jmean[i*L+j] 

    for i in range(L):
        hi_ = np.sum(h)
        h_IG[i,:] = h_IG[i,:] - hi_ + CorrectJ[i,:]

    return h_IG

def outputestimator(h_data,h,J_data,J):
    #------ Coupling Constant ------# 
    fname_data = "J_data_L"+str(L)+"_IsingGauge.dat"
    fname_model = "J_model_estimator_L"+str(L)+"_IsingGauge.dat"
    #fname = "J_model_estimator_CD1_L"+str(L)+"_IsingGauge.dat"
    
    f_data = open(fname_data, "w")
    f_model = open(fname_model, "w")
    #J_data_vec = np.reshape( np.copy(J_data), q**2*L**2)
    #J_model_vec = np.reshape( np.copy(J_model), q**2*L**2)
    for i in range(L):
        for j in range(L):
            if(j!=i):
                for a in range(q):
                    for b in range(q):
                        f_data.write(str(J_data[i*L+j][a*q+b])+"\n")
                        f_model.write(str(J_model[i*L+j][a*q+b])+"\n")
    f_data.close();f_model.close()

    #------ Local Field ------# 
    fname_data_h = "h_data_L"+str(L)+"_IsingGauge.dat"
    fname_model_h = "h_model_estimator_L"+str(L)+"_IsingGauge.dat"
    #fname = "h_model_estimator_CD1_L"+str(L)+"_IsingGauge.dat"
    
    f_data_h = open(fname_data_h, "w")
    f_model_h = open(fname_model_h, "w")
    h_data_vec = np.reshape( np.copy(h_data), q*L)
    h_model_vec = np.reshape( np.copy(h_model), q*L)
    
    for l in range(len(h_data_vec)):
        f_data_h.write(str(h_data_vec[l])+"\n")
    f_data_h.close()
    
    for l in range(len(h_model_vec)):
        f_model_h.write(str(h_model_vec[l])+"\n")
    f_model_h.close()

if __name__ == "__main__":
    fJname_data = "J_data_L"+str(L)+".dat"
    fJname_model = "J_model_estimator_L"+str(L)+".dat"
    #fJname_model = "J_model_estimator_CD1_L"+str(L)+".dat"
    
    fhname_data = "h_data_L"+str(L)+".dat"
    fhname_model = "h_model_estimator_L"+str(L)+".dat"
    
    J_data = J_init(fJname_data)
    J_model= J_init(fJname_model)
    h_data = h_init(fhname_data)
    h_model= h_init(fhname_model)
  
    J_IG_data = Ising_Gauge_J(J_data)
    h_IG_data = Ising_Gauge_h(h_data)
    
    J_IG_model = Ising_Gauge_J(J_model)
    h_IG_model = Ising_Gauge_h(h_model)
    outputestimator(h_IG_data, h_IG_model, J_IG_data, J_IG_model)


