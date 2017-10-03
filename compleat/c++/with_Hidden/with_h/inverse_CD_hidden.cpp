// Maximamu Likelihood estimation for inverse potts model. 
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <random>
using namespace std;

int q=3, L=4;
int p=L;
double Xi0 = 0.1, h0 = 0.1;

/*****************************/
vector<int> X(L,0);
vector<double> H(p,0);
vector < vector<double> > Xi( L*p, vector<double>(q, Xi0) );  
//vector < vector<double> > J_data( L*L, vector<double>(q*q,0) );  
vector < vector<double> > h( L, vector<double>(q,h0) );  
vector < vector<double> > h_data( L, vector<double>(q,0) );  
vector < vector<double> > Psi( L*p, vector<double>(q,0) );  
vector < vector<double> > Psi_model( L*p, vector<double>(q,0) );  
vector < vector<double> > Psi_h( L, vector<double>(q,0.0) );  
vector < vector<double> > Psi_model_h( L, vector<double>(q,0) );  

/*****************************/
int T_epoch = 2e1;
int T_interval = 1e2;
int N_sample = 1e5;
double lr = 1.0, lr_h = 1.0;
double eps = 1e-10; 
double error[2] ={0,0};
int K = 1;
/*****************************/

inline double std_rand()
{
    return rand() / (RAND_MAX + 1.0);
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

void init_Psi(){
	for(int i=0;i<L;++i){
	for(int a=0;a<q;++a){
		Psi_h[i][a] = 0.0;
		for(int m=0;m<p;++m){
			Psi[i*p+m][a] = 0.0;	
		}	
	}}
}

void init_Psi_model(){
	for(int i=0;i<L;++i){
	for(int a=0;a<q;++a){
		Psi_model_h[i][a] = 0.0;
		for(int m=0;m<p;++m){
			Psi_model[i*p+m][a] = 0.0;	
		}	
	}}
}

bool average_data_Hidden(){
	init_Psi();
	//string fname_in = "Psia_data_L"+to_string(L)+"_c++.dat";
	string fname_in = "test_training_data_L"+to_string(L)+"_c++.dat";
  	ifstream file(fname_in);
  	string line;
	vector<double> x_0(p,0);
	bool flag = true;
	int count = 0; 
	int a,b;
	vector<int> A(L,0);
	if (file.is_open())
	{
		count = 0; 
		while (getline (file,line,' ') && count<(N_sample*L) ) 
		{
			//X[count%L] = stoi(line);
			A[count%L] = stoi(line);
			if(count%L==L-1){
				fill(x_0.begin(), x_0.end(), 0);		
				//rather x_0 thatn H is proper.	
				/*x_0 correspond to  hidden variables*/	
				for(int m=0;m<p;++m){
					x_0[m] = 0.0; // NOTE: Necessary!
					for(int i=0;i<L;++i){
						a = int(A[i]);
						x_0[m] += Xi[i*p+m][a]; 
					}
					x_0[m] /= double(L);	
				}

				for(int i=0;i<L;++i){
					a = int(A[i]);
					Psi_h[i][a] += 1.0;
					//cout << "Psi_h["<<i<<"]["<<a<<"]="<<Psi_h[i][a] << endl;
					for(int m=0;m<p;++m){
						Psi[i*p+m][a] += x_0[m];
						//Psi[i*p+m][a] += H[m];
				}}
			}
			count  += 1;
			} //end while
		file.close();
	}//end if 
  	else{ cout << "Unable to open file" << endl;  
  	cout << fname_in << endl;
    	flag = false;
  	}
		
	for(int i=0;i<L;++i){
	for(int a=0;a<q;++a){
		Psi_h[i][a] /= double(N_sample);	
		//cout << "Norm: Psi_h["<<i<<"]["<<a<<"]="<<Psi_h[i][a] << endl;
		for(int m=0;m<p;++m){
			Psi[i*p+m][a] /= double(N_sample);	
	}}}	
	
	return flag;
}

bool average_model_CD_Hidden(){
	init_Psi_model();
	string fname_in = "test_training_data_L"+to_string(L)+"_c++.dat" ;
  	ifstream file(fname_in);
  	string line;
	vector<double> x_0(p,0);
	bool flag = true;
	int count = 0; 
	int a,b;
	if (file.is_open())
	{
		count = 0; 
		while (getline (file,line,' ') && count<(N_sample*L) ) 
		{
			X[count%L] = stoi(line);
			if(count%L==L-1){
				
				for(int k=0;k<K;++k){
					sampling_hidden();  //update H
					sampling_visible(); //update X	
				}
				
				for(int  i=0;i<L;++i){
					a = int(X[i]);
					Psi_model_h[i][a] += 1.0;
					for(int m=0;m<p;++m){
						Psi_model[i*p+m][a] += H[m];	
				}}
			}
			count  += 1;
			} //end while
		file.close();
	}//end if 
  	else{ cout << "Unable to open file" << endl;  
  	cout << fname_in << endl;
    	flag = false;
  	}
	
	for(int i=0;i<L;++i){
	for(int a=0;a<q;++a){
		Psi_model_h[i][a] /= double(N_sample);	
		for(int m=0;m<p;++m){
			Psi_model[i*p+m][a] /= double(N_sample);	
	}}}	
	
	return flag;
}

void gradient_descent_Hidden(){
	double dPsi,dPsi_h;
	int i,m,a;
	error[0]=0.0; error[1]=0.0;
	for(i=0;i<L;++i){
	for(a=0;a<q;++a){
            dPsi_h = (Psi_h[i][a]-Psi_model_h[i][a]);
            error[1] += dPsi_h*dPsi_h;
            h[i][a] += lr_h * (dPsi_h - eps*h[i][a]);
	    //cout << "Psi_h[i][a]=" << Psi_h[i][a]<< ", Psi_model_h[i][a]=" << Psi_model_h[i][a]<< endl;
		
		for(m=0;m<p;++m){
                	dPsi = (Psi[i*p+m][a] -  Psi_model[i*p+m][a]);
			//cout << "dPsi = " << dPsi << endl;
			error[0] += dPsi*dPsi;
                	Xi[i*p+m][a] += lr * (dPsi - eps*Xi[i*p+m][a]) ;
			//cout << "Psi[i*p+m][a]=" << Psi[i*p+m][a] << ", Psi_model[i*p+m][a]=" <<Psi_model[i*p+m][a] << endl;
		}
	}}
    	error[0] = sqrt( error[0] / float(L*p*q) )  ;
    	error[1]= sqrt( error[1] / float(L*q) )  ;
}

void outputestimator(){
    	// Coupling Constant 
	string fname = "Xi_model_estimator_CD_L"+to_string(L)+"_Hidden_c++.dat";
	ofstream file(fname);
	string fname_h = "h_model_estimator_CD_L"+to_string(L)+"_Hidden_c++.dat";
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
	
void output_statisticalmodel(){
	/****** 2nd cumulant *****/
	string fname_fij = "fij_CD_L"+to_string(L)+"_Hidden_c++.dat";
	string fname_pij = "pij_CD_L"+to_string(L)+"_Hidden_c++.dat";
	ofstream file_fij(fname_fij);
	ofstream file_pij(fname_pij);
	
	/******  1nd cumulant *****/
	string fname_fi = "fi_CD_L"+to_string(L)+"_Hidden_c++.dat";
	string fname_pi = "pi_CD_L"+to_string(L)+"_Hidden_c++.dat";
	ofstream file_fi(fname_fi);
	ofstream file_pi(fname_pi);
	
	for(int i=0;i<L;i++){
	for(int a=0;a<q;++a){
		file_fi << Psi_h[i][a] << endl;
		file_pi << Psi_model_h[i][a] << endl;
		for(int m=0;m<p;++m){
			file_fij << Psi[i*p+m][a] << endl; 
			file_pij << Psi_model[i*p+m][a] << endl; 
		}}}
	file_fij.close(); file_pij.close();
}

int main(){
	double time_start,time_end,dtime;
	int t;
	time_start = time(0);

	for(t=0;t<T_epoch;++t){
		average_data_Hidden();
		average_model_CD_Hidden();
		gradient_descent_Hidden();
		cout << t <<"  " << error[0] <<"  " << error[1] << endl;
	}
	time_end = time(0); 
	dtime = time_end - time_start;
	cout << "#ML: computational time = " << dtime <<endl;
	outputestimator();
	output_statisticalmodel();
}
