// Generative model of potts. 
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
using namespace std;

int q=3, L=4;
double J0 = 1.0;

/*****************************/
vector<int> X(L,0);
vector < vector<double> > J( L*L, vector<double>(q*q,0) );  
vector < vector<double> > h( L, vector<double>(q,0) );  
/*****************************/
int T_equil = 1e3;
int T_interval = 3e2;
int N_sampling = 1e5;
/*****************************/
int K = 1;
/*****************************/

inline double std_rand()
{
    return rand() / (RAND_MAX + 1.0);
}

bool J_init(){
	bool flag=true;
	//string fname_in = "J_model_estimator_L"+to_string(L)+"_c++.dat";
	string fname_in = "J_model_estimator_CD_L"+to_string(L)+"_c++.dat";
	//string fname_in = "J_data_CD_L"+to_string(L)+"_c++.dat";
  	ifstream file(fname_in);
  	string line;
  	int count = 0; int i=0,j=0,a,b,var; // this initialize is necessary.
	if (file.is_open()){
    	while ( getline (file,line,'\n'))
    	{
		/* NOTE
		 This method doesn't take into acount a case in which p = q
		 * */
		var = count;
	    	b = var % q; 
		var = var - b;
		a = (var/q)%q;
		var = (var - a*q)/(q*q);
		j = var%L;	
		i = int(var/L);
		J[L*i+j][a*q+b] = stod(line);
		count += 1;
    	}
    	file.close();
	flag = true;
	}
  	else{ cout << "Unable to open file" << endl;  
  	cout << fname_in << endl;
    	flag = false;
  	}
	return flag;
}

bool h_init(){
	bool flag=true;
	//string fname_in = "h_data_L"+to_string(L)+"_c++.dat";
	//string fname_in = "h_model_estimator_L"+to_string(L)+"_c++.dat";
	string fname_in = "h_model_estimator_CD_L"+to_string(L)+"_c++.dat";
	//string fname_in = "h_data_CD_L"+to_string(L)+"_c++.dat";
  	ifstream file(fname_in);
  	string line;
  	int count = 0; int i=0,a,var; // this initialize is necessary.
  	if (file.is_open()){
    	while ( getline (file,line,'\n'))
    	{
		var = count;
	    	a = var % q;
		var = var - a;
		i = int(var/q);	
		h[i][a] = stod(line);
		count += 1;
    	}
    	file.close();
	flag = true;
	}
  	else{ cout << "Unable to open file" << endl;  
  	cout << fname_in << endl;
    	flag = false;
  	}
	return flag;
}

double E_local(int i){
	int a = X[i],b;
	double E_i = 0.0;
	for(int j=0;j<L;++j){
		b = X[j];
		E_i += -J[i*L+j][a*q+b];	
	}
	E_i += -h[i][a];
	return E_i;	
}

double E(){
	double E_tot = 0.0;
	for(int i=0;i<L;++i){
		E_tot = E_local(i);
	}
	return E_tot / 2.0;
}

bool Metropolis(int i){
	double E_i, E_i_trial, dE, w;
	int Xi_save;
	bool accepted;
	
	E_i= E_local(i);
	Xi_save = X[i];
	X[i] = (X[i]+int(std_rand()*q+1)) % q;
	E_i_trial = E_local(i);
	dE = E_i_trial - E_i;
	w =exp(-dE);
       	accepted = false;
	if(w > std_rand())
       		accepted = true;
	else
		X[i] = Xi_save;
	return accepted;
}

int MonteCarlo_sweep(){
	double accepted;
	int i,l;
	for(l=0;l<L;++l){
		i = int(L*std_rand());
		if(Metropolis(i))
			accepted += 1;
	}
	return accepted;
}

int main(){
	J_init(); h_init();
	//string fname = "genmodel_data_ML_L"+to_string(L)+"_c++.dat" ;
	string fname = "genmodel_data_CD_L"+to_string(L)+"_c++.dat" ;
	ofstream file(fname);
	int t;
	t = 0;
	while(t<T_equil){
		MonteCarlo_sweep();	
		t+=1;
	}
	t = 0;
	while(t<T_interval*N_sampling){
		MonteCarlo_sweep();
		if(t%T_interval == 0){
		
		for(int i=0;i<L;++i){
			file << X[i] << " ";	
		}
			file << endl;
	}
		t+=1;
	}
	file.close();
}
