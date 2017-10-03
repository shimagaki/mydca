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
/*****************************/
int T_equil = 1e3;
int T_interval = 1e2;
int N_sampling = 1e4;
/*****************************/
inline double std_rand()
{
    return rand() / (RAND_MAX + 1.0);
}

void J_generalized_potts(){
	double r;
	for(int i=0;i<L;i++){
	for(int j=i+1;j<L;j++){
	//for(int j=0;j<L;j++){
		for(int a=0;a<q;++a){
		for(int b=0;b<q;++b){
			r = std_rand();
		//	J[i*L+j][a*q+b] = count;
			J[i*L+j][a*q+b] = r;
			J[j*L+i][b*q+a] = r;
		}}}}
	string fname = "J_data_L"+to_string(L)+"_c++.dat";
	ofstream file(fname);
	for(int i=0;i<L;i++){
	for(int j=0;j<L;j++){
		for(int a=0;a<q;++a){
		for(int b=0;b<q;++b){
			file << J[i*L+j][a*q+b] << endl; 
		}}}}
	file.close();
}
double E_local(int i){
	int a = X[i],b;
	double E_i = 0.0;
	for(int j=0;j<L;++j){
		b = X[j];
		E_i += -J[i*L+j][a*q+b];	
	}
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
	J_generalized_potts();
	string fname = "test_training_data_L" + to_string(L) + "_c++.dat";
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


