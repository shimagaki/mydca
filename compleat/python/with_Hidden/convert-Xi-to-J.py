# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt

q, L = 3, 4 
p = L
X = np.random.choice(q,L) # Chose L-variables from q-state .
Xi = np.zeros((L*p,q))
J = np.zeros((L**2, q**2))

def Xi_set():
    global Xi
    fname_in = "Xi_model_estimator_CD_L"+str(L)+"_Hidden.dat"
    fin = open(fname_in, "r")
     
    a=0; m=0; i=0; 
    for line in fin:
        line = line.replace('\n','') # do not need to use 
        item = line.split(' ')
        
        Xi[i*p+m][a]= float(item[0])
        a += 1
        if(a==q):
            a = 0 
            m += 1
            if(m==p):
                m = 0
                i += 1
    fin.close()

def convert_Xi_to_J():
    global J, Xi
    
    fname = "J_of_Xi_CD_L"+str(L)+".dat"
    fout = open(fname,"w")
    
    for i in range(L):
        for j in range(L):
            for a in range(q):
                for b in range(q):
                    for m in range(p):
                        if(i !=j ):
                            J[i*L+j][a*q+b] += Xi[i*p+m][a]*Xi[j*p+m][b]
                    fout.write(str(J[i*L+j][a*q+b])+"\n") 
    fout.close()

if __name__ == "__main__":
    Xi_set()
    convert_Xi_to_J()

