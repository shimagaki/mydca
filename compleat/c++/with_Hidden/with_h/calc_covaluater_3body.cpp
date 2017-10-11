#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
using namespace std;

int q=21, L=23;
int LL=L*L, qq=q*q;
//NOTE: N_sample is should be changed if you change input data X.
int N_sample_data = 340711; 
int N_sample_model =1e5; 
/*****************************/
vector<int> X(L,0);
vector < vector<double> > C_data( L*L, vector<double>(q*q,0) );  
vector < vector<double> > C3_data( L*L*L, vector<double>(q*q*q,0) );  
vector < vector<double> > f1_data( L, vector<double>(q,0) );  
vector < vector<double> > f2_data( L*L, vector<double>(q*q,0) );  
vector < vector<double> > f3_data( L*L*L, vector<double>(q*q*q,0) );  

vector < vector<double> > C_model( L*L, vector<double>(q*q,0) );  
vector < vector<double> > C3_model( L*L*L, vector<double>(q*q*q,0) );  
vector < vector<double> > f1_model( L, vector<double>(q,0) );  
vector < vector<double> > f2_model( L*L, vector<double>(q*q,0) );  
vector < vector<double> > f3_model( L*L*L, vector<double>(q*q*q,0) );  


/*****************************/

bool calc_C(){
	bool flag_data=true;
	string fname_in_data = "PF00096.dat";
  	ifstream file_in_data(fname_in_data);
  	string line_data;
  	int count_data = 0; 
	int a,b,c;//int i=0,j=0,k=0,a,b,c;
	double scale_ij,scale_ijk;
	if (file_in_data.is_open())
	{
		count_data = 0;
		while (getline (file_in_data,line_data,' ') && count_data<(N_sample_data*L) ) 
		{
			X[count_data%L] = stoi(line_data);
			if(count_data%L==0 && count_data > 0){
				for(int i=0;i<L;++i){
					a = X[i];
				       	f1_data[i][a] += 1.0/N_sample_data;
				for(int j=0;j<L;++j){ // Take into account diagonal parts.
					
					b = X[j];
					if(j==i){
						scale_ij=1.0/2.0; // 2!=2	
					}else{
						scale_ij=1.0;
					}
					f2_data[i*L+j][a*q+b] += scale_ij/N_sample_data;
				
					for(int k=0;k<L;++k){
					
						c = X[k];
						if(i==j && j==k){
							scale_ijk=1.0/6.0; // 3!=6
						}
						else if((i==j && j !=k) || (j==k && k!=i) || (k==i && i!=j)){
							scale_ijk=1.0/2.0; // 2!=2	
						}else{
							scale_ijk=1.0;
						}
						f3_data[i*LL+j*L+k][a*qq+b*q+c] += scale_ijk/N_sample_data ;
					}// end k-loop
				}// end j-loop
				}}
			count_data += 1;
		}
		cout << "count_data = " <<  count_data << endl;
		file_in_data.close();
		flag_data = true;
	}
  	else{ cout << "Unable to open file" << endl;  
  	cout << fname_in_data << endl;
    	flag_data = false;
  	}
	
	bool flag_model=true;
	string fname_in_model = "test_training_data_L"+to_string(L)+"_Hidden_c++.dat";
  	ifstream file_in_model(fname_in_model);
  	string line_model;
  	int count_model = 0; 
	//scale_ij,scale_ijk;
	
	if (file_in_model.is_open())
	{
		count_model = 0;
		while (getline (file_in_model,line_model,' ') && count_model<(N_sample_model*L) ) {
			X[count_model%L] = stoi(line_model);
			if(count_model%L==0 && count_model > 0){
				for(int i=0;i<L;++i){
					a = X[i];
				       	f1_model[i][a] += 1.0/N_sample_model;
				for(int j=0;j<L;++j){ // Take into account diagonal parts.
					
					b = X[j];
					if(j==i){
						scale_ij=1.0/2.0; // 2!=2	
					}else{
						scale_ij=1.0;
					}
					f2_model[i*L+j][a*q+b] += scale_ij/N_sample_model;
					
					for(int k=0;k<L;++k){
						
						c = X[k];
						if(i==j && j==k){
							scale_ijk=1.0/6.0; // 3!=6
						}
						else if((i==j && j !=k) || (j==k && k!=i) || (k==i && i!=j)){
							scale_ijk=1.0/2.0; // 2!=2	
						}else{
							scale_ijk=1.0;
						}
						f3_model[i*LL+j*L+k][a*qq+b*q+c] += scale_ijk/N_sample_model ;

					}// end k-loop	
				}//end j-loop
				}}
			count_model += 1;
		}
 		cout << "count_model = " <<  count_model << endl;
		file_in_model.close();
		flag_model = true;
	}
  	else{ cout << "Unable to open file" << endl;  
  	cout << fname_in_model << endl;
    	flag_model = false;
  	}
	bool flag=false;	
	if(flag_data && flag_model){
		flag = true;
	}
	
	for(int i=0;i<L;++i){
	for(int j=0;j<L;++j){
		for(int a=0;a<q;++a){
		for(int b=0;b<q;++b){
			C_data[i*L+j][a*q+b] = f2_data[i*L+j][a*q+b] - f1_data[i][a]*f1_data[j][b];
			
			C_model[i*L+j][a*q+b] = f2_model[i*L+j][a*q+b] - f1_model[i][a]*f1_model[j][b];
			
			for(int k=0;k<L;++k){
				for(int c=0;c<q;++c){
					
					C3_data[i*LL+j*L+k][a*qq+b*q+c] = f3_data[i*LL+j*L+k][a*qq+b*q+c]
							-(f1_data[i][a]*f2_data[j*L+k][b*q+c]
							+ f1_data[j][b]*f2_data[k*L+i][c*q+a]
							+ f1_data[k][c]*f2_data[i*L+j][a*q+b] )
							+2*f1_data[i][a]*f1_data[j][b]*f1_data[k][c];
					
					C3_model[i*LL+j*L+k][a*qq+b*q+c] = f3_model[i*LL+j*L+k][a*qq+b*q+c]
							-(f1_model[i][a]*f2_model[j*L+k][b*q+c]
							+ f1_model[j][b]*f2_model[k*L+i][c*q+a]
							+ f1_model[k][c]*f2_model[i*L+j][a*q+b] )
							+2*f1_model[i][a]*f1_model[j][b]*f1_model[k][c];
				}
			}
		}}}}
		
	string fname_out_data = "Correlation_Teacher_L"+to_string(L)+"3_c++.dat";
	ofstream file_out_data(fname_out_data);
	string fname_out_model = "Correlation_CD_L"+to_string(L)+"_c++.dat";
	ofstream file_out_model(fname_out_model);
	
	string fname_out3_data = "3Correlation_Teacher_L"+to_string(L)+"_c++.dat";
	ofstream file_out3_data(fname_out3_data);
	string fname_out3_model = "3Correlation_CD_L"+to_string(L)+"_c++.dat";
	ofstream file_out3_model(fname_out3_model);
	
	for(int i=0;i<L;++i){
	for(int j=0;j<L;++j){
		
		for(int a=0;a<q;++a){
		for(int b=0;b<q;++b){
			file_out_data << C_data[i*L+j][a*q+b] << endl; 
			file_out_model << C_model[i*L+j][a*q+b] << endl; 
		}}
		
		for(int k=0;k<L;++k){
			for(int c=0;c<q;++c){
				file_out3_data << C3_data[i*LL+j*L+k][a*qq+b*q+c] << endl; 
				file_out3_model << C3_model[i*LL+j*L+k][a*qq+b*q+c]  << endl; 
		}}
	}}
       	file_out_data.close();
       	file_out_model.close();
	
	return flag;
}

int main(){
	calc_C();
}

