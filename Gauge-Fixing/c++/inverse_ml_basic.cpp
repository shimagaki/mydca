// Maximamu Likelihood estimation for inverse potts model. 
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
using namespace std;

int q=3, L=5;
double J0 = 0.5, h0 = 0.3;

/*****************************/
vector<int> X(L,0);
vector < vector<double> > J( L*L, vector<double>(q*q,J0) );  
vector < vector<double> > J_data( L*L, vector<double>(q*q,0) );  
vector < vector<double> > h( L, vector<double>(q,h0) );  
vector < vector<double> > h_data( L, vector<double>(q,0) );  
vector < vector<double> > Psi( L*L, vector<double>(q*q,0) );  
vector < vector<double> > Psi_model( L*L, vector<double>(q*q,0) );  
vector < vector<double> > Psi_h( L, vector<double>(q,0) );  
vector < vector<double> > Psi_model_h( L, vector<double>(q,0) );  

/*****************************/
int T_epoch = 1e3;
int T_equil = 1e3;
int T_interval = 1e2;
int N_sample = 1e4;
int M = 1e2 ; // number of sample of average by statistical model. 
double lr = 0.001, lr_h = 0.001;
double eps = 0.0;
double error[2] ={0,0};
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
	E_i += -h[i][a];
	return E_i;	
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
	//string fname_in = "J_data_L"+to_string(L)+".dat";
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
bool h_data_init(){
	bool flag=true;
	string fname_in = "h_data_L"+to_string(L)+"_c++.dat";
	//string fname_in = "h_data_L"+to_string(L)+".dat";
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
		h_data[i][a] = stod(line);
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

void J_model_init(){
	for(int i=0;i<L;++i){
		for(int a=0;a<q;++a){
		for(int b=0;b<q;++b){
                	J[i*L+i][a*q+b] = 0; 
	}}}
}

void h_model_init(){
	for(int i=0;i<L;++i){
		for(int a=0;a<q;++a){
                	h[i][a] = 0; 
	}}
}

bool average_data_J(){
	bool flag=true;
	string fname_in = "test_training_data_L"+to_string(L)+"_c++.dat";
	//string fname_in = "test_training_data_L"+to_string(L)+".dat";
  	ifstream file(fname_in);
  	string line;
  	int count = 0; 
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
					Psi[i*L+j][a*q+b] += 1.0/N_sample;
					Psi[j*L+i][b*q+a] += 1.0/N_sample;
		}}}
			count += 1;
		}
		cout << "count = " <<  count << endl;
		file.close();
		flag = true;
	}
  	else{ cout << "Unable to open file" << endl;  
  	cout << fname_in << endl;
    	flag = false;
  	}
	return flag;
}
bool average_data_h(){
	bool flag=true;
	string fname_in = "test_training_data_L"+to_string(L)+"_c++.dat";
	//string fname_in = "test_training_data_L"+to_string(L)+".dat";
  	ifstream file(fname_in);
  	string line;
  	int count = 0; 
	int i,a;
	if (file.is_open())
	{
		count = 0;
		while (getline (file,line,' ') && count<(N_sample*L) ) 
		{
			X[count%L] = stoi(line);
			if(count%L==0 && count > 0){
				for(i=0;i<L;++i){
					a = X[i];
					Psi_h[i][a] += 1.0/N_sample;
				}}
			count += 1;
		}
		file.close();
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
void Psi_model_h_init(){
	for(int i=0;i<L;++i){
		for(int a=0;a<q;++a){
			Psi_model_h[i][a] = 0.0;
	}}	
}

void X_init(){
	for(int i=0;i<L;++i){
		X[i] = int(q*std_rand());
	}
}

int  average_model(){
	Psi_model_init();
	Psi_model_h_init();
	X_init();
	int t,i,j,a,b,accepted=0;
	for(t=0;t<T_equil;++t){
		MonteCarlo_sweep();
	}

	for(t=0;t<T_interval*M;++t){
		accepted += MonteCarlo_sweep();
		if(t%T_interval == 0){
			for(i=0;i<L;++i){
				a = X[i];
			       	Psi_model_h[i][a] += 1.0/M;
				for(j=i+1;j<L;++j){
					b = X[j];
					Psi_model[i*L+j][a*q+b] += 1.0/M; 
					Psi_model[j*L+i][b*q+a] += 1.0/M;
	}}}}
	return accepted ;
}

void gradient_descent(){
	double dPsi,dPsi_h;
	int i,j,a,b;
	error[0]=0.0; error[1]=0.0;
	for(i=0;i<L;++i){
	for(j=i+1;j<L;++j){
		for(a=0;a<q;++a){
		for(b=0;b<q;++b){
			
			dPsi = (Psi[i*L+j][a*q+b] -  Psi_model[i*L+j][a*q+b]);
			//cout << "dPsi = " << dPsi << endl;
			error[0] += dPsi*dPsi;
			J[i*L+j][a*q+b] +=  lr * (dPsi - eps*J[i*L+j][a*q+b]); 
			
			dPsi = (Psi[j*L+i][b*q+a] - Psi_model[j*L+i][b*q+a]);
			error[0] += dPsi*dPsi ;
			J[j*L+i][b*q+a] +=  lr * (dPsi - eps*J[j*L+i][b*q+a]); 
		}}	
	}}
    	error[0] = sqrt(error[0]) / (L*L*q*q);
	for(i=0;i<L;++i){
		for(a=0;a<q;++a){
			dPsi_h = (Psi_h[i][a]-Psi_model_h[i][a]);
			h[i][a] += lr_h * (dPsi_h - eps*h[i][a]);
			error[1] += dPsi_h*dPsi_h;
	}}
    	error[1]= sqrt(error[1]) / (L*q);

}

void outputestimator(){
    	// Coupling Constant 
	string fname = "J_model_estimator_L"+to_string(L)+"_c++.dat";
	ofstream file(fname);
	for(int i=0;i<L;i++){
	for(int j=0;j<L;j++){
		for(int a=0;a<q;++a){
		for(int b=0;b<q;++b){
			file << J[i*L+j][a*q+b] << endl; 
		}}}}
	file.close();
	
    	// Magnetic Field 
	string fname_h = "h_model_estimator_L"+to_string(L)+"_c++.dat";
	ofstream file_h(fname_h);
	for(int i=0;i<L;i++){
		for(int a=0;a<q;++a){
			file_h << h[i][a] << endl; 
		}}
	file_h.close();
}

int main(){
	double time_start,time_end,dtime;
	int t;
	J_data_init();
	J_model_init();
	average_data_J();
	average_data_h();
	time_start = time(0);
	for(t=0;t<T_epoch;++t){
		average_model();
		gradient_descent();
		cout << t <<"  " << error[0] <<"  " << error[1] << endl;
	}
	time_end = time(0); 
	dtime = time_end - time_start;
	cout << "#ML: computational time = " << dtime <<endl;
	outputestimator();
	
}
