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
int N_sample = 1e3;
/*****************************/
vector<int> X(L,0);
vector < vector<double> > C( L*L, vector<double>(q*q,0) );  
vector < vector<double> > f2( L*L, vector<double>(q*q,0) );  
vector < vector<double> > f1( L, vector<double>(q,0) );  

/*****************************/

bool calc_C(){
	bool flag=true;
	//string fname_in = "test_training_data_L"+to_string(L)+"_c++.dat";
	string fname_in = "genmodel_data_ML_L"+to_string(L)+"_c++.dat";
  	ifstream file_in(fname_in);
  	string line;
  	int count = 0; 
	int i=0,j=0,a,b;
	if (file_in.is_open())
	{
		count = 0;
		while (getline (file_in,line,' ') && count<(N_sample*L) ) 
		{
			X[count%L] = stoi(line);
			if(count%L==0 && count > 0){
				for(int i=0;i<L;++i){
					a = X[i];
				       	f1[i][a] += 1.0/N_sample;
				//for(int j=i+1;j<L;++j){
				for(int j=i;j<L;++j){ // Take into account diagonal parts.
					b = X[j];
					f2[i*L+j][a*q+b] += 1.0/N_sample;
					f2[j*L+i][b*q+a] += 1.0/N_sample;
		}}}
			count += 1;
		}
		cout << "count = " <<  count << endl;
		file_in.close();
		flag = true;
	}
  	else{ cout << "Unable to open file" << endl;  
  	cout << fname_in << endl;
    	flag = false;
  	}
	
	for(int i=0;i<L;++i){
	for(int j=i+1;j<L;++j){
		for(int a=0;a<q;++a){
		for(int b=0;b<q;++b){
			C[i*L+j][a*q+b] = f2[i*L+j][a*q+b] - f1[i][a]*f1[j][b];
			C[i*L+j][a*q+b] = C[i*L+j][a*q+b];
		}}}} 
		
	//string fname_out = "Correlation_Teacher_L"+to_string(L)+"_c++.dat";
	string fname_out = "Correlation_ML_L"+to_string(L)+"_c++.dat";
	ofstream file_out(fname_out);
	
	for(int i=0;i<L;++i){
	for(int j=0;j<L;++j){
		for(int a=0;a<q;++a){
		for(int b=0;b<q;++b){
			file_out << C[i*L+j][a*q+b] << endl; 
		}}}}
       	file_out.close();	
	
	return flag;
}

int main(){
	calc_C();
}




