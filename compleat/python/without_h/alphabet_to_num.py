# -*- coding: utf-8 -*- 
import numpy as np
import matplotlib.pyplot as plt
L = 23
code_alpha = ["A","C","D","E","F","G","H","I",
"K","L","M","N","P","Q","R","S","T","V","W","Y","-"]

def read_write():
    #fname_in = "test_data.dat"
    fname_in = "PF00096_msa.txt"
    fname_out = "test_output.dat"

    fin = open(fname_in, "r")
    fout = open(fname_out, "w")
    for line in fin:
        line = line.replace('\n','')
        if(line[0]!=">"): 
            item = line.split()
            for l in range(L):
                num = code_alpha.index(line[l]) 
                fout.write(str(num)+" ")  
            fout.write("\n") 
    fin.close()
    fout.close()

if __name__ == "__main__":
    read_write()















