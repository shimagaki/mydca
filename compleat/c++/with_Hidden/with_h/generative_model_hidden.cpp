// Generative model of potts. 
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <random>
#include <numeric>
using namespace std;

int q=3, L=16;
int p=L;
double Xi0=0.1, h0 = 1.0;

/*****************************/
vector<int> X(L,0);
vector<double> H(p,0);
vector < vector<double> > Xi( L*p, vector<double>(q, Xi0) );  
vector < vector<double> > h( L, vector<double>(q,0) );  
//vector < vector<double> > Psi( L*p, vector<double>(q,0) );  
//vector < vector<double> > Psi_h( L, vector<double>(q,0.0) );  
/*****************************/
int T_equil = 1e3;
int T_interval = 1e2;
//int N_sampling = 1e5;
int N_sampling = 1e4;
/*****************************/
int K = 1;
/*****************************/

inline double std_rand()
{
    return rand() / (RAND_MAX + 1.0);
}

void X_set(){
	for(int i=0;i<L;++i){
		X[i] = int(q*std_rand()) ;
	}
}

void set_h_Xi(){
	for(int i=0;i<L;++i){
		for(int a=0;a<q;++a){
			h[i][a] = h0 * std_rand();
			for(int m=0;m<p;++m){
				Xi[i*p+m][a] = Xi0 * std_rand();
			
	}}}
}

bool Xi_init(){
	bool flag=true;
	string fname_in = "Xi_model_estimator_CD_L"+to_string(L)+"_Hidden_c++.dat";
  	ifstream file(fname_in);
  	string line;
  	int i=0,m=0,a=0; 
	if (file.is_open()){
    	while ( getline (file,line,'\n'))
    	{
		Xi[i*p+m][a] = stod(line);
		m += 1;
		if(m==p){//NOTE: m moves faster than a.
			a+=1; m=0;	
			if(a==q){
				i+=1; a=0;
			}
		}		
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
	string fname_in = "h_model_estimator_CD_L"+to_string(L)+"_Hidden_c++.dat";
  	ifstream file(fname_in);
  	string line;
  	int count = 0; int i=0,a=0,var; // this initialize is necessary.
		
  	if (file.is_open()){
    	while ( getline (file,line,'\n'))
    	{
		h[i][a] = stod(line);
		a += 1;
		if(a==q){
			i+=1; a=0;
		}	
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

void get_normaldist(vector<double> mean, double std){
	random_device rd;
	mt19937 gen(rd());
	for(int m=0;m<p;++m){
		normal_distribution<> d(mean[m],std);
		H[m] = d(gen);
		//cout << m << " " << mean[m] << " " << H[m] << endl;	
	}		
}

int my_sampling_gibbsdist(vector<double> prob){
	random_device rd;
	mt19937 gen(rd());
	using dist_type = discrete_distribution<>;
	discrete_distribution<> d(prob.begin(), prob.end());
	int resolute = d(gen);
	return resolute; 
}

double E_local_hidden(int i,int a){
	double E_i = 0.0;
	for(int m=0;m<p;++m){
		E_i += - Xi[i*p+m][a]*H[m];
	}
	E_i += -h[i][a];
	return E_i;	
}

void sampling_hidden(){
	int a;
	vector<double> x_0(p,0);
	for(int m=0;m<p;++m){
		for(int i=0;i<L;++i){
			a = int(X[i]);
			x_0[m] += Xi[i*p+m][a];
		}
		x_0[m] /= double(L);
	}
	get_normaldist(x_0, 1.0/L);// mean=x_0, std=1/L, size=p
}

void sampling_visible(){
	double E_min,pdis_sum;
	vector<double> pdis(q,0);
	vector<double> E_vec(q,0);
	
	for(int i=0;i<L;++i){
	 	for(int a=0;a<q;++a){	E_vec[a] = E_local_hidden(i,a);	}
		E_min = *min_element(E_vec.begin(), E_vec.end());
		pdis_sum = 0.0;
		for(int a=0;a<q;++a){
			pdis[a] = exp( - (E_local_hidden(i,a) - E_min) );	
			pdis_sum += pdis[a];
		}
		//pdis_sum = accumulate(pdis.begin(), pdis.end(), 0);
	 	for(int a=0;a<q;++a){ pdis[a] /= pdis_sum ; }
		X[i] = my_sampling_gibbsdist(pdis);	
	}
}

void outputestimator(){
    	// Coupling Constant 
	string fname = "Xi_training_L"+to_string(L)+"_Hidden_c++.dat";
	ofstream file(fname);
	string fname_h = "h_training_L"+to_string(L)+"_Hidden_c++.dat";
	ofstream file_h(fname_h);
	
	for(int i=0;i<L;i++){
	for(int a=0;a<q;++a){
		file_h << h[i][a] << endl; 
		for(int m=0;m<p;++m){
			file << Xi[i*p+m][a] << endl; 
		}}}
	file.close(); 
	file_h.close();
}

int main(){
	h_init();
	Xi_init();	
	
	/* 
	//set_h_Xi();
	//outputestimator();	
	*/
	
	//string fname = "genmodel_data_ML_L"+to_string(L)+"_c++.dat" ;
    	string fname  = "test_training_data_L"+to_string(L)+"_Hidden_c++.dat";
	ofstream file(fname);
	int t;
	t = 0;
	//X_set();
	while(t<T_equil){
		sampling_hidden();
		sampling_visible();
		t += 1;	
	}	
	
	t = 0;
	while(t<N_sampling*T_interval){
		sampling_hidden();
		sampling_visible();
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
