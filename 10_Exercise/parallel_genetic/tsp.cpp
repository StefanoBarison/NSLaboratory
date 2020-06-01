#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>

//Useful libraries
#include<vector>
#include<algorithm>
#include<iterator>
#include<string>
#include "mpi.h"

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

	int x=(int)rnd->Rannyu(1,size-n);
	int y=(int)rnd->Rannyu(1,size-n);

	while(abs(x-y)<n){
		y=(int)rnd->Rannyu(1,size-n);
	}

	for(int i=0;i<n;i++){
		swap(_chromosome[x+i],_chromosome[y+i]);
	}
}

void individual:: Uniform_swap_mutate(double p_u, Random* rnd){
	int size=_chromosome.size();

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

void individual:: Inversion_mutate(int n,Random * rnd){
	int size=_chromosome.size();

	int x= (int)rnd->Rannyu(1,size-n);
	int y=x+n;

	reverse(_chromosome.begin()+x,_chromosome.begin()+y);

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

		for(int j=0;j<1000;j++){
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

void population:: Wheel_selection(map cities, Random* rnd){
		// First, order the cities
	this-> Sort(cities);

	//Then, generate the vector of probabilities, the probabilities are generated proportionally to the highet possible path
	vector<double> fitness;

	for(int i=0;i<_size;i++){
		fitness.push_back(_pop[i].Get_lenght());
	}

	//double tot=accumulate(fitness.begin(),fitness.end(),0.00);
	double max=_pop[_size-1].Get_lenght();

	for(int i=0;i<_size;i++){
		//fitness[i]=1.0-fitness[i]/tot;
		fitness[i]=1.0-fitness[i]/max;	
	}

	vector<double> probabilities;
	for(int i=0;i<_size;i++){
		double t=accumulate(fitness.begin(),fitness.begin()+i,0.00);
		probabilities.push_back(t);
	}

	for(int i=0;i<_size;i++){
		probabilities[i]/=probabilities[_size-1];
		//cout<<probabilities[i]<<",";
	}

	//cout<<endl;
	//Now let's perform the biased roulette wheel selection

	vector<individual> new_pop;
	

	for(int j=0;j<_size;j++){
		double r=rnd->Rannyu();

		for(int k=0;k<_size-1;k++){
			if(r> probabilities[k] && r<probabilities[k+1] ){
				new_pop.push_back(_pop[k]);
			}
		}
	}

	//Substitute the old generation with the new
	_pop=new_pop;

	//Evaluate the new population

	for(int i=0;i<_size;i++){
		_pop[i].Evaluate(cities);
	}
}

void population::Elitism(map cities, int elite, Random* rnd){

	// Order the cities, and save the best candidates
	this-> Sort(cities);
	int n_elite=elite;

	vector<individual> new_pop;
	for(int i=0;i<n_elite;i++){
		new_pop.push_back(_pop[i]);
	}

	//Then, generate the vector of probabilities, the probabilities are generated with respect to the longest possible path
	vector<double> fitness;

	for(int i=0;i<_size;i++){
		fitness.push_back(_pop[i].Get_lenght());
	}

	//double tot=accumulate(fitness.begin(),fitness.end(),0.00);
	double max=_pop[_size-1].Get_lenght();

	for(int i=0;i<_size;i++){
		//fitness[i]=1.0-fitness[i]/tot;
		fitness[i]=1.0-fitness[i]/max;	
	}

	vector<double> probabilities;
	for(int i=0;i<_size;i++){
		double t=accumulate(fitness.begin(),fitness.begin()+i,0.00);
		probabilities.push_back(t);
	}

	for(int i=0;i<_size;i++){
		probabilities[i]/=probabilities[_size-1];
		//cout<<probabilities[i]<<",";
	}

	//cout<<endl;
	//Now let's perform the biased roulette wheel selection

	for(int j=0;j<_size-n_elite;j++){
		double r=rnd->Rannyu();
		//cout<<r<<",";

		for(int k=0;k<_size-1;k++){
			if(r> probabilities[k] && r<probabilities[k+1] ){
				new_pop.push_back(_pop[k]);
			}
		}
	}

	//Substitute the old generation with the new
	_pop=new_pop;

	//Evaluate the new population

	for(int i=0;i<_size;i++){
		_pop[i].Evaluate(cities);
	}
}

void population:: Evolutive_step(map cities,string type,Random* rnd){
	
	if(type=="Roulette"){
		//cout<<"Performing roulette wheel selection..."<<endl;
		//Mutate the individuals with a probability m_p and crossover with a probability c_p
		double m_p1=0.03;
		double m_p2=0.09;
		double m_p3=0.10;
		double c_p=0.60;
		
		for(int i=0;i<_size-1;i++){
			double r1=rnd->Rannyu();
			if(r1<c_p){
				int r_i=(int)rnd->Rannyu(0,_size);
				Crossover(&_pop[i],&_pop[r_i],rnd);
			}
		}


		for(int i=0;i<_size;i++){
			double r2=rnd->Rannyu();
			if(r2<m_p1){
				this->Get_individual(i)->Swap_mutate(rnd);
			}
			
			if(r2>m_p1 && r2<m_p2){
				this->Get_individual(i)->Push_back_mutate(2);
			}
			
			if(r2>m_p2 && r2<m_p3){
				this->Get_individual(i)->Multi_swap_mutate(3,rnd);
			}
			/*if(r2>m_p3 && r2<m_p4){
				this->Get_individual(i)->Inversion_mutate(4,rnd);
			}*/
		}

		//Now select the best individuals and take to the next generation

		this->Wheel_selection(cities,rnd);
	}


	else if(type=="Elitism"){
		//cout<<"Performing roulette wheel selection with elitism..."<<endl;
		//First, we have to save the best candidates

		vector<individual> elite;
		int n_elite=5;

		this->Sort(cities);

		for(int i=0;i<n_elite;i++){
			elite.push_back(_pop[i]);
		}

		//Now perform mutations and crossovers
		double m_p1=0.05;
		double m_p2=0.10;
		double m_p3=0.12;
		double m_p4=0.15;
		double c_p=0.50;
		
		for(int i=0;i<_size-1;i++){
			double r1=rnd->Rannyu();
			if(r1<c_p){
				int r_i=(int)rnd->Rannyu(0,_size);
				Crossover(&_pop[i],&_pop[r_i],rnd);
			}
		}


		for(int i=0;i<_size;i++){
			double r2=rnd->Rannyu();
			if(r2<m_p1){
				this->Get_individual(i)->Swap_mutate(rnd);
			}
			if(r2>m_p1 && r2<m_p2){
				this->Get_individual(i)->Push_back_mutate(2);
			}

			if(r2>m_p2 && r2<m_p3){
				this->Get_individual(i)->Multi_swap_mutate(3,rnd);
			}

			if(r2>m_p3 && r2<m_p4){
				this->Get_individual(i)->Inversion_mutate(5,rnd);
			}
		}

		this->Sort(cities);
		//Now re introduce the elite eliminating the longest paths
		for(int i=0;i<n_elite;i++){
			_pop[_size-1-i]=elite[i];
		}

		this->Sort(cities);

		//...and perform selection
		this->Elitism(cities,n_elite,rnd);
	}
}

void population:: Simulated_annealing(map cities,double temp, Random* rnd){
	this -> Sort(cities);
	// In simulated annealing we modify the solution candidates with random mutation (maybe uniform mutation)
	// in order to obtain a new candidate, we evaluate its path lenght and decide to accept or reject it

	//usually simulated annealing is performed on only one solution, but here we can do for _size individuals
cout<<"=================="<<endl;
	for(int i=0;i<_size;i++){
		cout<<"Individual: "<<i<<endl;
		_pop[i].Print();
		_pop[i].Print_lenght();

		vector<int> proposal=_pop[i].Get_genes();

		//Create a new individual

		individual new_ind(proposal);
		individual old_ind(proposal);

		//Modify randomly

		double m_p1=0.30, m_p2=0.40, m_p3=0.60, m_p4=0.95, m_p5=1.00;
		double p_u=0.01;
		double c_p=0.6;
		double r_p=rnd->Rannyu();
		double r_c=rnd->Rannyu();

		cout<<"Number for random mutation: "<<r_p<<endl;

		if(_pop[i].Get_lenght()>60.0){
		if(r_p<m_p1) new_ind.Swap_mutate(rnd);
		else if(r_p<m_p2 && r_p>=m_p1) new_ind.Push_back_mutate(3);
		else if(r_p<m_p3 && r_p>=m_p2) new_ind.Multi_swap_mutate(2,rnd);
		else if(r_p<m_p4 && r_p>=m_p3) new_ind.Inversion_mutate(10,rnd);
		else if(r_p<m_p5 && r_p>=m_p4) new_ind.Uniform_swap_mutate(p_u,rnd);
		}
		
		//Crossover between the old and the new proposal
		if(r_c<c_p){
			Crossover(& old_ind,& new_ind,rnd);
			cout<<"Crossover"<<endl;
		}

		new_ind.Evaluate(cities);
		old_ind.Evaluate(cities);

		if(old_ind<new_ind){
			new_ind=old_ind;
		}

		cout<<"New proposal: "<<endl;
		new_ind.Print();
		new_ind.Print_lenght();

		//Now accept/reject the move

		if(new_ind<_pop[i]){
			_pop[i]=new_ind;
			cout<<"***Accepted***"<<endl;
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
			cout<<"Random number: "<<b_r<<endl;

			if(b_r<prob){
				_pop[i]=new_ind;
			}

			else if(b_r>prob){
				cout<<"***rejected***"<<endl;
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


void Crossover(individual * t1,individual * t2, Random* rnd){
	vector<int> p1=t1->Get_genes();
	vector<int> p2=t2->Get_genes();

	int size=p1.size();

	vector<int> s1;
	vector<int> s2;


	int r=(int)rnd->Rannyu(1,size);

	//cout<<r<<endl;

	for(int i=0;i<r;i++){
		s1.push_back(p1[i]);
		s2.push_back(p2[i]);
	}

	for(int i=0;i<size;i++){

		int a=p1[i];
		int b=p2[i];

		if(any_of(p1.begin()+r,p1.end(),[b](int k){return k==b;})){
			s1.push_back(p2[i]);
		}
		if(any_of(p2.begin()+r,p2.end(),[a](int k){return k==a;})){
			s2.push_back(p1[i]);
		}
	}


	t1->Set_genes(s1);
	t2->Set_genes(s2);
}

