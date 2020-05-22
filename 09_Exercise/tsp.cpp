#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>

//Useful libraries
#include<vector>
#include<algorithm>
#include<iterator>
#include<random>

//Personal
#include "tsp.h"


// Classe individual

using namespace std;

void individual:: Evaluate(map cities){
	double sum=0.0;
	int size= _chromosome.size();
	for(int i=0;i<size-1;i++){
		sum+=cities.Distance(_chromosome[i],_chromosome[i+1]);
	}

	sum+=cities.Distance(_chromosome[size-1],0);

	_lenght=sum;
}

double individual::Get_lenght(){
	return _lenght;
}

void individual::Print(){

	int size =_chromosome.size();
	for(int i=0;i<size;i++){
		cout<<_chromosome[i]<<" ";
	}
	cout<<endl;
}

void individual::Print_lenght(){
	cout<<"Tour's lenght: "<<_lenght<<endl;
}

vector<int> individual::Get_genes(){
	return _chromosome;
}

void individual::Set_genes(vector<int> new_genes){
	_chromosome=new_genes;
}

bool individual:: Check(){
	vector<int> reference;
	int size=_chromosome.size();

	for(int i=0;i<size;i++){
		reference.push_back(i);
	}

	bool check=is_permutation(_chromosome.begin(),_chromosome.end(),reference.begin());

	return check;
}

individual::~individual(){
	_chromosome.clear();
}

// Classe population

void population:: Initialize(){
	vector<int> reference;
	for(int i=0;i<_ncities;i++){
		reference.push_back(i);
	}

	for(int i=0;i<_size;i++){
		random_shuffle(reference.begin()+1,reference.end());
		individual new_one(reference);

		_pop.push_back(new_one);
	}

}

individual* population::Get_individual(int i){
	return &_pop[i];
}


double population::Get_size(){
	return _size;
}





// Classe city

double city::Get_x(){
	return _x;
}

double city::Get_y(){
	return _y;
}

void city::Print(){
	cout<<_x<<","<<_y<<endl;
}

city:: ~city(){}
// Classe map





double map::Distance(int i,int j){
	return _distances[i][j];
}

double map:: Get_size(){
	return _ncities;
}

void map:: Create_d_matrix(){
	for(int i=0;i<_ncities;i++){
		vector<double> row;
		for(int j=0;j<_ncities;j++){
			double d= sqrt(pow(_cities[i].Get_x()-_cities[j].Get_x(),2)+pow(_cities[i].Get_y()-_cities[j].Get_y(),2));
			row.push_back(d);

		}
		_distances.push_back(row);
	}
}


void map::Circle_initialize(){
	double l=10.0;
	random_device radm;
    mt19937_64 mt(radm());
    uniform_real_distribution<> distribution(0, 2*M_PI);

    for(int i=0;i<_ncities;i++){
    	double theta=distribution(mt);

    	city new_city(l*cos(theta),l*sin(theta));

    	_cities.push_back(new_city);
    }

}

void map::Print(){
	for(int i=0;i<_ncities;i++){
		_cities[i].Print();
	}
}


void map::Print_d_matrix(){
	int wd=10;

	for(int i=0;i<_ncities;i++){
		for(int j=0;j<_ncities;j++){
			cout<<setw(wd)<<_distances[i][j];
		}
		cout<<endl;
	}
}

