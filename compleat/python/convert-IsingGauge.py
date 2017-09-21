# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
L, q = 2, 3 

def J_init(fJname):
    
    fin = open(fJname, "r")
    J_vec = [] 
    i = 0
    for line in fin:
        line = line.replace('\n','') 
        J_vec.append(float(line)) 
        i += 1
    ###NOTE#######################NOTE###
    J = np.reshape( np.copy(J_vec), (L**2,q**2) )
    fin.close()
    return J

def h_init(fhname):
    fin = open(fhname, "r")
    h_vec = [] 
    i = 0
    for line in fin:
        line = line.replace('\n','')
        h_vec.append(float(line)) 
        i += 1
    h = np.reshape( np.copy(h_vec), (L,q) )
    fin.close()
    return h

def Ising_Gauge_J(J):
    J_IG = np.copy(J)
    for i in range(L):
        for j in range(1+i,L):
            Jij_b = np.zeros(q) # sum of col
            Jija_ = np.zeros(q) # sum of row 
                 
            for c in range(q):
                for d in range(q):
                    Jij_b[c] += J[i*L+j][d*q+c]/q
                    Jija_[c] += J[i*L+j][c*q+d]/q
                    
            Jij__ = sum(Jij_b) / q 
            
            for a in range(q):
                for b in range(q):
                    J_IG[i*L+j][a*q+b] = J[i*L+j][a*q+b] \
                            - Jij_b[b] - Jija_[a] + Jij__
    
                    J_IG[j*L+i][b*q+a] = J_IG[i*L+j][a*q+b]
    return J_IG

def Ising_Gauge_h(h,J):
    h_IG = np.copy(h)
    CorrectJ = np.zeros((L**2,q)) 
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
        
    for i in range(L):
        hi_ = np.sum(h)
        h_IG[i,:] = h_IG[i,:] - hi_ + CorrectJ[i,:]

    return h_IG

def outputestimator(J_data,J):
    #------ Coupling Constant ------# 
    fname_data = "J_data_L"+str(L)+"_IsingGauge.dat"
    fname = "J_model_estimator_L"+str(L)+"_IsingGauge.dat"
    #fname = "J_model_estimator_CD1_L"+str(L)+"_IsingGauge.dat"
    f_data = open(fname_data, "w")
    f = open(fname, "w")
    J_data_vec = np.reshape( np.copy(J_data), q**2*L**2)
    J_model_vec = np.reshape( np.copy(J), q**2*L**2)
    
    for l in range(len(J_data_vec)):
        f_data.write(str(J_data_vec[l])+" \n")
    f_data.close()
    
    for l in range(len(J_model_vec)):
        f.write(str(J_model_vec[l])+" \n")
    f.close()

if __name__ == "__main__":
    fJname_data = "J_data_L"+str(L)+".dat"
    fJname_model = "J_model_estimator_L"+str(L)+".dat"
    #fJname_model = "J_model_estimator_CD1_L"+str(L)+".dat"
    
    J_data = J_init(fJname_data)
    J_model= J_init(fJname_model)
  
    J_IG = Ising_Gauge_J(J_data)
    #h_IG = Ising_Gauge_h(h_data,J_data)
    print "J_IG(J_data)= ", J_IG 
    
    print "\n"
    J_IG = Ising_Gauge_J(J_model)
    print "J_IG(J_model)= ", J_IG 
    outputestimator(J_data,J_model) 
