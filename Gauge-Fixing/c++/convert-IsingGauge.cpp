#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
using namespace std;
int q=3, L=8;
double J0=1.0, h0=1.0;

vector < vector<double> > J( L*L, vector<double>(q*q,J0) );  
vector < vector<double> > h( L, vector<double>(q,h0) );  
vector < vector<double> > J_Iging( L*L, vector<double>(q*q,J0) );  
vector < vector<double> > h_Iging( L, vector<double>(q,h0) );  
bool J_init(){
	bool flag=true;
	string fname_in = "J_data_L"+to_string(L)+"_c++.dat";
	//string fname_in = "J_model_estimator_L"+to_string(L)+"_c++.dat";
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
	string fname_in = "h_data_L"+to_string(L)+"_c++.dat";
	//string fname_in = "h_model_estimator_L"+to_string(L)+"_c++.dat";
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


void convert_J_Ising(){
        
	for(int i=0;i<L;i++){
	for(int j=i+1;j<L;j++){
		
		vector<double> Jij_b(q,0);  
		vector<double> Jija_(q,0);  
		double sum_b = 0.0, Jij__ =0.0, r, temp;
		
		for(int c=0;c<q;++q){
		for(int d=0;d<q;++q){
                   	temp =  J[i*L+j][d*q+c]/q;
			temp += Jij__/q;
			Jij_b[c] += temp;
                    	Jija_[c] += J[i*L+j][c*q+d]/q;
		}}
		
		for(int a=0;a<q;++a){	
		for(int b=0;b<q;++b){
                    	J_Iging[i*L+j][a*q+b] = J[i*L+j][a*q+b] - Jij_b[b] - Jija_[a] + Jij__;
                    	J_Iging[j*L+i][b*q+a] = J_Iging[i*L+j][a*q+b];
		}}
    
		}}
	string fname = "J_data_L"+to_string(L)+"_Ising_c++.dat";
	//string fname = "J_data_ML_L"+to_string(L)+"_Ising_c++.dat";
	//string fname = "J_data_CD_L"+to_string(L)+"_Ising_c++.dat";
	ofstream file(fname);
	for(int i=0;i<L;i++){
	for(int j=0;j<L;j++){
		for(int a=0;a<q;++a){
		for(int b=0;b<q;++b){
			file << J_Iging[i*L+j][a*q+b] << endl; 
		}}}}
	file.close();
}

void convert_h_Ising(){
	for(int i=0;i<L;++i){
	for(int j=0;j<L;++j){
		if(i !=j){
			for(int a=0;a<q;++a){
			for(int b=0;b<q;++b){
                        	CorrectJ[i*L+j][a] += J[i*L+j][a*L+b] / q 
			
			
			
			}}	
		
		}	
	
	}}	



}














	
int main(){
	J_init();
	h_init();
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
