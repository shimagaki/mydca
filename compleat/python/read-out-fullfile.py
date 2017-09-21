#!/usr/bin/env python                                                           
#coding:utf-8                                                                   
L = 2
name= "J_data_L"+str(L)+".dat"
f = open(name,'r')
Allf = f.read()

text = Allf.replace('\n','')
text = text.replace('\r','')
print len(text) 
print text[0]
f.close()

def J_init(fJname):
    
    fin = open(fJname, "r")
    J_vec = [] 
    i = 0
    for line in fin:
        #item = line.split(' ')
        line.replace('\n','')
        print "item=", line 
        #J_vec.append(map(float, line)) 
        J_vec.append(float(line)) 
        i += 1
    ###NOTE#######################NOTE###
    fin.close()

if __name__ == "__main__":

    fname = "J_data_L"+str(L)+".dat"
    J_init(fname)

 
