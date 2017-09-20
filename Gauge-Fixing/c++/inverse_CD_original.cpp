// Maximamu Likelihood estimation for inverse potts model. 
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
using namespace std;

int q=2, L=16;
double J0 = 1.0;

/*****************************/
vector<int> X(L,0);
vector < vector<double> > J( L*L, vector<double>(q*q,0) );  
vector < vector<double> > J_data( L*L, vector<double>(q*q,0) );  
vector < vector<double> > Psi( L*L, vector<double>(q*q,0) );  
vector < vector<double> > Psi_model( L*L, vector<double>(q*q,0) );  

/*****************************/
int T_epoch = 100;
int T_equil = 100;
int T_interval = 1e2;
///int N_sample = 1e4;
int N_sample = 1800;
int M = 100 ; // number of sample of average by statistical model. 
double lr = 0.05;
double eps = 0.0;
/*****************************/
int K = 1;
/*****************************/
inline double std_rand()
{
    return rand() / (RAND_MAX + 1.0);
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

bool J_data_init(){
	bool flag=true;
	string fname_in = "J_data_L"+to_string(L)+"_c++.dat";
  	ifstream file(fname_in);
  	string line;
  	int count = 0; int i=0,j=0,a,b,var; // this initialize is necessary.
  	if (file.is_open()){
    	while ( getline (file,line,'\n'))
    	{
		var = count;
	    	b = var % q;
		var = var - b;
		a = (var/q)%q;
		var = (var - a*q)/(q*q);
		j = var%L;	
		i = int(var/L);
		J_data[L*i+j][a*q+b] = stod(line);
		count += 1;
		/*
		cout << J_data[L*i+j][a*q+b] << " ";
		if(count%L == 0)
			cout << endl;
		*/
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

bool average_data(){
	bool flag=true;
	string fname_in = "test_training_data_L"+to_string(L)+"_c++.dat";
  	ifstream file(fname_in);
  	string line;
  	int count; 
	int i=0,j=0,a,b;
	if (file.is_open())
	{
		count = 0;
		while (getline (file,line,' ') && count<(N_sample*L) ) 
		{
			X[count%L] = stoi(line);
			if(count%L==0 && count > 0){
				for(int i=0;i<L;++i){
				for(int j=i+1;j<L;++j){
					a = X[i]; b = X[j];
					Psi[i*L+j][a*q+b] += 1.0;
					Psi[j*L+i][b*q+a] += 1.0;
		}}}
			count += 1;
		}
		for(i=0;i<L;++i){
		for(j=i+1;j<L;++j){
			for(int a=0;a<q;++a){
			for(int b=0;b<q;++b){
				Psi[i*L+j][a*q+b] /= N_sample;
				Psi[j*L+i][b*q+a] /= N_sample;
			}}}}
		file.close();
		flag = true;
	}
  	else{ cout << "Unable to open file" << endl;  
  	cout << fname_in << endl;
    	flag = false;
  	}
	return flag;
}

void Psi_model_init(){
	for(int i=0;i<L;++i){
	for(int j=0;j<L;++j){
		for(int a=0;a<q;++a){
		for(int b=0;b<q;++b){
			Psi_model[i*L+j][a*q+b] = 0.0;
			Psi_model[j*L+i][b*q+a] = 0.0;
	}}}}	
}

void X_init(){
	for(int i=0;i<L;++i){
		X[i] = int(q*std_rand());
	}
}

int  average_model_CD(){
	Psi_model_init();
	X_init();
	int t,i,j,a,b,accepted=0,count;
	bool flag = true;
	
	string fname_in = "test_training_data_L"+to_string(L)+"_c++.dat";
  	ifstream file(fname_in);
  	string line;
	if (file.is_open())
	{
		count = 0;
		while (getline (file,line,' ') && count<(N_sample*L) ) 
		{
			X[count%L] = stoi(line);
			if(count%L==0 && count > 0 ){
				/*contrastive divergence.*/
            			for(int k=0;k<K;k++){
					MonteCarlo_sweep();
				}
				
				/*taking an average.*/
				for(int i=0;i<L;++i){
				for(int j=i+1;j<L;++j){
					a = X[i]; b = X[j];
					Psi_model[i*L+j][a*q+b] += 1.0/N_sample;
					Psi_model[j*L+i][b*q+a] += 1.0/N_sample;
				}}
			}	
			count += 1;
		}
		file.close();
		flag = true;
	}
  	else{ cout << "Unable to open file" << endl;  
  	cout << fname_in << endl;
    	flag = false;
  	}
	return accepted ;
}

double gradient_descent(){
	double error = 0.0,dPsi;
	int i,j,a,b;
	for(i=0;i<L;++i){
	for(j=i+1;j<L;++j){
		for(a=0;a<q;++a){
		for(b=0;b<q;++b){
			//cout << Psi[i*L+j][a*q+b] << "  " << Psi_model[i*L+j][a*q+b] << endl;
			dPsi = (Psi[i*L+j][a*q+b] -  Psi_model[i*L+j][a*q+b]);
			error += dPsi*dPsi ;
			J[i*L+j][a*q+b] +=  lr * (dPsi - eps*J[i*L+j][a*q+b]); 
			
			dPsi = (Psi[j*L+i][b*q+a] - Psi_model[j*L+i][b*q+a]);
			error += dPsi*dPsi ;
			J[j*L+i][b*q+a] +=  lr * (dPsi - eps*J[j*L+i][b*q+a]); 
		}}	
	}}
    error = sqrt(error) / (L*L*q*q);
	return error;
}

void outputestimator(){
	//string fname = "J_model_estimator_L"+to_string(L)+"_c++.dat";
	string fname = "J_model_estimator_L"+to_string(L)+"_c++.dat";
	ofstream file(fname);
	for(int i=0;i<L;i++){
	for(int j=0;j<L;j++){
		for(int a=0;a<q;++a){
		for(int b=0;b<q;++b){
			file << J[i*L+j][a*q+b] << endl; 
		}}}}
	file.close();
}

int main(){
	double error,time_start,time_end,dtime;
	int t;
	bool flag;

	flag = J_data_init();
	average_data();
	time_start = time(0);
	for(t=0;t<T_epoch;++t){
		average_model_CD();
		error = gradient_descent();
		cout << t <<"  " << error << endl;
	}
	time_end = time(0); 
	dtime = time_end - time_start;
	cout << "#ML: computational time = " << dtime <<endl;
	outputestimator();
	
	
}
