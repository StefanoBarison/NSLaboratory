#include "tsp.h"

#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<random>


using namespace std;


int main (){

	map a(10);

	a.Circle_initialize();

	a.Print();

	a.Create_d_matrix();

	a.Print_d_matrix();

	population p(10,a);


	p.Initialize();

	for(int i=0;i<10;i++){
		p.Get_individual(i)->Evaluate(a);
		p.Get_individual(i)->Print();
		p.Get_individual(i)->Print_lenght();
		
	}


	return 0;
	
}