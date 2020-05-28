#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>

//Useful libraries
#include<vector>
#include<algorithm>
#include<iterator>
#include<string>

//Personal
#include "tsp.h"
#include "Random.h"

   


using namespace std;

/////////////////////
// Class individual//
/////////////////////


bool individual:: operator ==(individual &ind2){
if(this->Get_lenght() == ind2.Get_lenght()){
		return true;
	}

	else{
		return false;
	}
	
}
bool individual::operator <(individual & ind2){
	if(this->Get_lenght() < ind2.Get_lenght()){
		return true;
	}

	else{
		return false;
	}
}

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
		cout<<setw(3)<<_chromosome[i];
	}
	cout<<endl;
}

void individual::Print_lenght(){
	cout<<_lenght<<endl;
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

void individual::Swap_mutate(Random* rnd){
	int size=_chromosome.size();
	int x=(int)rnd->Rannyu(1,size);
	int y=(int)rnd->Rannyu(1,size);

	while(x==y){
		y=(int)rnd->Rannyu(1,size);
	}

	swap(_chromosome[x],_chromosome[y]);
}

void individual:: Push_back_mutate(int n){
	vector<int> new_chromo;

	new_chromo.push_back(0);

	for(unsigned int i=n+1;i<_chromosome.size();i++){
		new_chromo.push_back(_chromosome[i]);
	}

	for(int i=1;i<=n;i++){
		new_chromo.push_back(_chromosome[i]);
	}

	_chromosome=new_chromo;
}

void individual:: Multi_swap_mutate(int n, Random* rnd){
	int size=_chromosome.size();

	int x=(int)rnd->Rannyu(1,size-n+1);
	int y=(int)rnd->Rannyu(1,size-n+1);

	while(abs(x-y)<n){
		y=(int)rnd->Rannyu(1,size-n+1);
	}

	for(int i=0;i<n;i++){
		swap(_chromosome[x+i],_chromosome[y+i]);
	}
}

void individual:: Uniform_swap_mutate(double p_u, Random* rnd){
	int size=_chromosome.size();
	uniform_int_distribution<> dis(1,size-1);
	uniform_real_distribution<> r_dist(0,1);

	for(int i=1;i<size;i++){
		double r=rnd->Rannyu();
		if(r<p_u){
			int y=(int)rnd->Rannyu(1,size);
			while(y==i){
				y=(int)rnd->Rannyu(1,size);
			}
			swap(_chromosome[i],_chromosome[y]);
		}
	}

}

////////////////////
//Class population//
////////////////////

void population:: Initialize(Random* rnd){
	vector<int> reference;
	for(int i=0;i<_ncities;i++){
		reference.push_back(i);
	}

	for(int i=0;i<_size;i++){
		individual new_one(reference);

		for(int j=0;j<100;j++){
			new_one.Swap_mutate(rnd);
			new_one.Uniform_swap_mutate(0.3,rnd);
		}

		_pop.push_back(new_one);
	}

}

void population:: Evaluate_all(map cities){
	for(int i=0;i<_size;i++){
		_pop[i].Evaluate(cities);
	}

}

individual* population::Get_individual(int i){
	return &_pop[i];
}


double population::Get_size(){
	return _size;
}

void population::Print_best(map cities){
	this-> Sort(cities);

	cout<<"Best candidate :"<<endl;
	_pop[0].Print();
	_pop[0].Print_lenght();
}


individual* population:: Get_best(map cities){

	this->Sort(cities);
	return &_pop[0];
}

double population:: Get_best_average(map cities){
	this->Sort(cities);
	int half=_size/2;

	double sum=0.0;

	for(int i=0;i<half;i++){
		sum+=_pop[i].Get_lenght();
	}

	sum/=half;

	return sum;
}

void population:: Sort(map cities){

	for(int i=0;i<_size;i++){
		_pop[i].Evaluate(cities);
	}

	sort(_pop.begin(),_pop.end());

}

void population:: Simulated_annealing(map cities,double temp, Random* rnd){
	this -> Sort(cities);
	// In simulated annealing we modify the solution candidates with random mutation (maybe uniform mutation)
	// in order to obtain a new candidate, we evaluate its path lenght and decide to accept or reject it

	//usually simulated annealing is performed on only one solution, but here we can do for _size individuals

	for(int i=0;i<_size;i++){
		vector<int> proposal=_pop[i].Get_genes();

		//Create a new individual

		individual new_ind(proposal);

		//Modify randomly and evaluate

		double m_p1=0.10, m_p2=0.50, m_p3=0.7, m_p4=1.00;
		double p_u=0.001;
		double r_p=rnd->Rannyu();

		
		if(r_p<m_p1) new_ind.Swap_mutate(rnd);
		else if(r_p<m_p2 && r_p>=m_p1) new_ind.Push_back_mutate(2);
		else if(r_p<m_p3 && r_p>=m_p2) new_ind.Multi_swap_mutate(6,rnd);
		else if(r_p<m_p4 && r_p>=m_p3) new_ind.Uniform_swap_mutate(p_u,rnd);

		new_ind.Evaluate(cities);

		//Now accept/reject the move

		if(new_ind<_pop[i]){
			_pop[i]=new_ind;
		}

		else if(new_ind==_pop[i]){
			_pop[i]=new_ind;
		}

		else if(_pop[i]<new_ind){
			double e1=new_ind.Get_lenght();
			double e2=_pop[i].Get_lenght();
			double prob= exp((e2-e1)/temp);
			cout<< "Boltzmann probability: "<<prob<<endl;

			double b_r=rnd->Rannyu();

			if(b_r<prob){
				_pop[i]=new_ind;
			}
		}
	}
}


//////////////
//Class city//
//////////////

double city::Get_x(){
	return _x;
}

double city::Get_y(){
	return _y;
}

void city::Print(){
	cout<<setw(10)<<_x<<setw(10)<<_y<<endl;
}

city:: ~city(){}


/////////////
//Class map//
/////////////

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


void map::Circle_initialize(Random* rnd){
	double l=10.0;

    for(int i=0;i<_ncities;i++){
    	double theta=rnd->Rannyu(0,2*M_PI);

    	city new_city(l*cos(theta),l*sin(theta));

    	_cities.push_back(new_city);
    }

}

void map::Square_initialize(Random* rnd){
	double l=10.0;


	for(int i=0;i<_ncities;i++){
		double x=rnd->Rannyu(-l,l);
		double y=rnd->Rannyu(-l,l);

		city new_city(x,y);
		_cities.push_back(new_city);
	}

}

void map::File_initialize(istream & indata){
	double x,y;

	while(!indata.eof()){
		indata>>x>>y;

		city new_city(x,y);
		_cities.push_back(new_city);
	}

	_ncities=_cities.size();
}

void map::Print_cities(){
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

