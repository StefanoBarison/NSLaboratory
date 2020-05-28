#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

//Personal
#include "tsp.h"
#include "Random.h"

using namespace std;

  
int main(int argc, char *argv[]){

  Random rnd;
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
    Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open()){
    while ( !input.eof() ){
      input >> property;
      if( property == "RANDOMSEED" ){
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed,p1,p2);
      }
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;

 //Now we have initialized the random number generator

 // We can start with our simulated annealing algorithm, define useful quantities
  int n_cities=32;
  int n_ind=1;
  int n_steps=150000;

  ifstream indata;
  indata.open("circle.dat");

  //initialize the map with the cities from previous lesson
  map m(n_cities);
  m.File_initialize(indata);
  m.Print_cities();
  m.Create_d_matrix();
  m.Print_d_matrix();

  //Create a population of a single individual to perform SA
  population p(n_ind,m);
  p.Initialize(& rnd);
  p.Evaluate_all(m);

  //Vector to save the results

  vector<double> best_result;

  //.. and a vector of temperatures

  vector<double> temperatures;

  for(int i=0;i<150000;i++){
    /*if(i<1000) temperatures.push_back(4.0);
    else if(i>=1000 && i<2000) temperatures.push_back(3.5);
    else if(i>=2000 && i<3000) temperatures.push_back(3.0);
    else if(i>=3000 && i<5000) temperatures.push_back(2.0);
    else if(i>=5000 && i<8000) temperatures.push_back(0.5);
    else if(i>=8000 && i<10000) temperatures.push_back(0.1);
    else if(i>=10000 && i<15000) temperatures.push_back(0.05);*/

    if(i<30000) temperatures.push_back(4.0-i*0.0001);
    else if(i>=30000 && i<50000) temperatures.push_back(0.03);
    else if(i>= 50000) temperatures.push_back(0.003);
    
  }

  //Perform the SA and print the results on a file

  ofstream outdata;
  outdata.open("test.dat");

  for (int i = 0; i < n_steps; ++i){

    cout<<"------------------------------------------------"<<endl;
    cout<<"Step "<<i+1<<" at temperature "<<temperatures[i]<<endl<<endl;
    p.Simulated_annealing(m,temperatures[i],& rnd);
    best_result.push_back(p.Get_individual(0)->Get_lenght());

    if(i==n_steps-1){
      p.Print_best(m);
    }
  }

  for(unsigned int i=0;i<best_result.size();i++){
    outdata<<i<<","<<best_result[i]<<endl;
  }


  return 0;
}
