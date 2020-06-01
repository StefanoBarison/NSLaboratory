#ifndef TSP
#define TSP

#include<iostream>
#include<fstream>

#define _USE_MATH_DEFINES
#include<cmath>

//Useful libraries

#include<vector>
#include<algorithm>
#include<iterator>
#include<random>
#include<string>
#include"mpi.h"
#include"Random.h"


class city{

public:
	city(double x, double y){
		_x=x;
		_y=y;
	}

	city():city(0,0){}

	~city();

	double Get_x();
	double Get_y();
	void Print();
private:
	double _x;
	double _y;
};


class map{

public:

	map(){
		_ncities=0;
	}

	map(int n){
		_ncities=n;
	}
	~map(){
		_cities.clear();
		_distances.clear();
	}

	void Circle_initialize(Random *rnd);
	void Square_initialize(Random *rnd);
	void File_initialize(std::istream & indata);
	void Create_d_matrix();

	double Distance(int i, int j);
	double Get_size();

	void Print_cities();

	void Print_d_matrix();


private:
	int _ncities;
	std::vector<std::vector<double>> _distances;
	std::vector<city> _cities;

};

class individual{
	
public:
	
	//individual();
	individual(std::vector<int> data){
		_chromosome=data;
		_lenght=0.0;
	}

	~individual();

	bool operator ==(individual & ind2);
	bool operator <(individual & ind2);

	void Evaluate(map cities);

	double Get_lenght();

	void Print();

	void Print_lenght();

	std::vector<int> Get_genes();

	void Set_genes(std::vector<int> new_genes);

	bool Check();

	void Swap_mutate(Random* rnd);

	void Push_back_mutate(int n);

	void Multi_swap_mutate(int n, Random *rnd);

	void Uniform_swap_mutate(double p_u, Random* rnd);

	void Inversion_mutate(int n, Random* rnd);

private:
	std::vector<int> _chromosome;
	double _lenght;

};


class population{

public:
	population(){
		_size=0;
		_ncities=0;
	}
	population(int n, map cities){
		_size=n;
		_ncities=cities.Get_size();
	}
	~population(){
		_pop.clear();
	}

	void Initialize(Random* rnd);
	void Evaluate_all(map cities);

	individual* Get_individual(int i);
	double Get_size(); 

	void Print_best(map cities);
	individual* Get_best(map cities);
	double Get_best_average(map cities);

	void Sort(map cities); //A function to sort the population from the shortest path to the longest

	void Wheel_selection(map cities, Random* rnd);
	void Elitism(map cities, int elite, Random* rnd);
	void Evolutive_step(map cities,std::string type,Random* rnd);

	void Simulated_annealing(map cities, double temp, Random* rnd); //A function to perform a simulated annealing search based 
													  // on a decreasing temperature and on the lenght difference
													  // between paths, it uses the Random number generator

private:
	int _size;
	int _ncities;
	std::vector<individual> _pop;

};


void Crossover(individual * t1,individual * t2, Random* rnd);

#endif