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

	void Circle_initialize();
	//void square_initialze();
	void Create_d_matrix();

	double Distance(int i, int j);
	double Get_size();

	void Print();

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
	}

	~individual();

	void Evaluate(map cities);

	double Get_lenght();

	void Print();

	void Print_lenght();

	std::vector<int> Get_genes();

	void Set_genes(std::vector<int> new_genes);

	bool Check();


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

	void Initialize();

	individual* Get_individual(int i);
	double Get_size();

private:
	int _size;
	int _ncities;
	std::vector<individual> _pop;

};



#endif