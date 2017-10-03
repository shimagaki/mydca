# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
L, q = 16, 3 
Jmean = []
Jmean_k = []

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
    global Jmean, Jmean_k
    J_IG = np.copy(J)
    
    Jmean = np.zeros(L**2)
    Jmean_k = np.zeros((L**2,q))
    for i in range(L):
        for j in range(L):
            Jij_b = np.zeros(q) # sum of col
            Jija_ = np.zeros(q) # sum of row 
                 
            for c in range(q):
                for d in range(q):
                    Jij_b[c] += J[i*L+j][d*q+c]/q
                    Jija_[c] += J[i*L+j][c*q+d]/q
                    Jmean[i*L+j] += J[i*L+j][d*q+c]/q**2
                Jmean_k[i*L+j][c] = Jija_[c]
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
    for i in range(L):
        for k in range(q):
            for j in range(L):
                if(j!=i):
                    CorrectJ[i][k] +=Jmean_k[i*L+j][k]-Jmean[i*L+j] 
                    #CorrectJ[i][k] +=Jmean_k[j*L+i][k]-Jmean[i*L+j] 

    for i in range(L):
        hi_ = np.sum(h[i,:])/q
        h_IG[i,:] = h_IG[i,:] - hi_ + CorrectJ[i,:]

    return h_IG

def outputestimator(h_data,h_model,J_data,J_model):
    #------ Coupling Constant ------# 
    fname_data = "J_data_L"+str(L)+"_IsingGauge_c++.dat"
    #fname_model = "J_model_estimator_L"+str(L)+"_IsingGauge_c++.dat"
    fname_model = "J_model_estimator_CD_L"+str(L)+"_IsingGauge_c++.dat"
    
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
    f_data.close()
    f_model.close()

    #------ Local Field ------# 
    fname_data_h = "h_data_L"+str(L)+"_IsingGauge_c++.dat"
    #fname_model_h = "h_model_estimator_L"+str(L)+"_IsingGauge_c++.dat"
    fname_model_h = "h_model_estimator_CD_L"+str(L)+"_IsingGauge_c++.dat"
    
    f_data_h = open(fname_data_h, "w")
    f_model_h = open(fname_model_h, "w")
    
    for i in range(L):
        for a in range(q):
            f_data_h.write(str(h_data[i][a])+"\n")
            f_model_h.write(str(h_model[i][a])+"\n")
    f_data_h.close()
    f_model_h.close()

if __name__ == "__main__":
    fJname_data = "J_data_L"+str(L)+"_c++.dat"
    #fJname_model = "J_model_estimator_L"+str(L)+"_c++.dat"
    fJname_model = "J_model_estimator_CD_L"+str(L)+"_c++.dat"
    
    fhname_data = "h_data_L"+str(L)+"_c++.dat"
    #fhname_model = "h_model_estimator_L"+str(L)+"_c++.dat"
    fhname_model = "h_model_estimator_CD_L"+str(L)+"_c++.dat"
    
    J_data = J_init(fJname_data)
    J_model= J_init(fJname_model)
    h_data = h_init(fhname_data)
    h_model= h_init(fhname_model)
  
    
    J_IG_data = Ising_Gauge_J(J_data)
    h_IG_data = Ising_Gauge_h(h_data)
    J_IG_model = Ising_Gauge_J(J_model)
    h_IG_model = Ising_Gauge_h(h_model)
    
    outputestimator(h_IG_model, h_IG_data, J_IG_model, J_IG_data)
    
    """ 
    sum_Jmod, sum_Jdat = 0.0, 0.0 
    print "#####"
    for i in range(L):
        for j in range(L):
            for c in range(q):
                Jij_b_dat, Jija__dat= 0.0, 0.0
                Jij_b_mod, Jija__mod= 0.0, 0.0
                for d in range(q):
                    sum_Jmod += J_IG_model[i*L+j][d*q+c] 
                    sum_Jdat += J_IG_data[i*L+j][d*q+c] 
                    Jij_b_dat += J_IG_data[i*L+j][d*q+c] 
                    Jija__dat += J_IG_data[i*L+j][c*q+d] 
                    Jij_b_mod += J_IG_model[i*L+j][d*q+c] 
                    Jija__mod += J_IG_model[i*L+j][c*q+d] 
                #print Jij_b_dat," ", Jija__dat, " ", Jij_b_mod, " ", Jija__mod
                print  Jij_b_mod, " ", Jija__mod
    print "sum_Jdat=", sum_Jdat,", sum_Jmod=", sum_Jmod
    print "#####"
    sum_hmod, sum_hdat = 0.0, 0.0
    for i in range(L):
        hi__dat = 0.0
        hi__mod = 0.0
        for b in range(q):
            sum_hmod+=h_IG_model[i][b]
            sum_hdat+=h_IG_data[i][b]
            hi__dat += h_IG_data[i][b]
            hi__mod += h_IG_model[i][b]
        print  hi__mod
    print "sum_hdat=", sum_hdat, ", sum_hmod=",sum_hmod
    """ 
