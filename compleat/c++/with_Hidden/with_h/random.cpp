#include <iostream>
#include <string>
#include <cmath>
#include <random>
#include <vector>
#include <map>
using namespace std;

int L = 16;
vector <double> r(L);
//map<int, int> hist;

void normal_rand(){
	random_device rd;
	mt19937 gen(rd());
	//normal_distribution<> d(0.0, 1.3);
	normal_distribution<> d(0.0,2.0);
	double r;
	for(int n=0; n<100; ++n) {
        	r = d(gen) ;	
		//++hist[ round(r) ];
		cout << r << endl;
	}
}

void a_random(){
	random_device rd;
	mt19937 gen(rd());
	using dist_type = discrete_distribution<>;
	
	vector<int> w(4);
       	w[0]=40,w[1]=10,w[2]=10;w[3]=40 ;	
	discrete_distribution<> d(w.begin(),w.end());	
    map<int, int> m;
    for(int n=0; n<10; ++n) {
        int resolute = d(gen);
	++m[resolute];
	cout << n << " " << resolute << endl;
    }
	int len_m = m.size();
	for(int i=0;i<len_m;++i){
		cout << i << " " << m[i] << endl;	
	}

}

int main(){
	//normal_rand();
	a_random();
	//int len_hist = hist.size();
	/*	
	for(int i=0;i<len_hist;++i){
		cout << i << " " << hist[i] << endl;	
	}
	*/	
	return 0;
}
